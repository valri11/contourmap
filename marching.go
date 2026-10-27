package contourmap

import (
	"math"
	"sync"
)

const closed = -math.MaxFloat64

// Edge slots within a cell.
const (
	eTop = iota
	eBottom
	eLeft
	eRight
)

const (
	noEdge   = -1 // no segment starts on this edge
	consumed = -2 // segment already emitted
)

// segment is a directed piece of contour inside one cell, from the edge it
// enters through to the edge it leaves through. Directions are chosen so the
// high side is always on the same hand, which makes every crossed edge the
// entry of exactly one segment and the exit of at most one.
type segment struct{ from, to uint8 }

var cellSegments = [16][]segment{
	1:  {{eTop, eLeft}},
	2:  {{eRight, eTop}},
	3:  {{eRight, eLeft}},
	4:  {{eLeft, eBottom}},
	5:  {{eTop, eBottom}},
	6:  {{eRight, eTop}, {eLeft, eBottom}},
	7:  {{eRight, eBottom}},
	8:  {{eBottom, eRight}},
	9:  {{eTop, eLeft}, {eBottom, eRight}},
	10: {{eBottom, eTop}},
	11: {{eBottom, eLeft}},
	12: {{eLeft, eRight}},
	13: {{eTop, eRight}},
	14: {{eLeft, eTop}},
}

// saddleHigh holds the alternative resolution of the two saddle cases, in
// which the high corners are connected through the cell.
var saddleHigh = [16][]segment{
	6: {{eLeft, eTop}, {eRight, eBottom}},
	9: {{eTop, eRight}, {eBottom, eLeft}},
}

func fraction(z0, z1, z float64) float64 {
	const eps = 1e-9
	var f float64
	if z0 == closed {
		f = 0
	} else if z1 == closed {
		f = 1
	} else if z0 != z1 {
		f = (z - z0) / (z1 - z0)
	}
	f = math.Max(f, eps)
	f = math.Min(f, 1-eps)
	return f
}

// cellPoint returns the crossing point on one edge of the cell whose
// upper-left sample is (x, y).
func cellPoint(slot uint8, x, y int, ul, ur, ll, lr, z float64) Point {
	fx, fy := float64(x), float64(y)
	switch slot {
	case eTop:
		return Point{fx + fraction(ul, ur, z), fy}
	case eBottom:
		return Point{fx + fraction(ll, lr, z), fy + 1}
	case eLeft:
		return Point{fx, fy + fraction(ul, ll, z)}
	default:
		return Point{fx + 1, fy + fraction(ur, lr, z)}
	}
}

func connectHigh(rule SaddleRule, ul, ur, ll, lr, z float64) bool {
	switch rule {
	case SaddleConnectHigh:
		return true
	case SaddleCenter, SaddleAsymptotic:
		if ul == closed || ur == closed || ll == closed || lr == closed {
			return false
		}
		if d := ul + lr - ur - ll; rule == SaddleAsymptotic && d != 0 {
			return (ul*lr-ur*ll)/d > z
		}
		return (ul+ur+ll+lr)/4 > z
	default:
		return false
	}
}

// scratch holds per-call working memory, indexed by edge. Horizontal edge
// (x, y)-(x+1, y) has index y*(w-1)+x; vertical edge (x, y)-(x, y+1) has index
// (w-1)*h + y*w + x. Between calls every next entry is noEdge and every hasIn
// entry is false, so a scratch can be reused without clearing it.
type scratch struct {
	next  []int32 // edge the segment entering through e leaves through
	hasIn []bool  // some segment leaves through e
	pts   []Point // crossing point on e; valid only for edges touched this call
	from  []int32 // entry edges in scan order
}

var scratchPool = sync.Pool{New: func() any { return new(scratch) }}

// maxPooledEdges bounds the scratch kept in the pool (about 44 MB, a grid of
// roughly 1000 x 1000). Larger scratches are left to the garbage collector so
// one big grid doesn't pin its buffers for every later small call.
const maxPooledEdges = 1 << 21

// maxEdges is the largest edge count an int32 edge id can address.
const maxEdges = math.MaxInt32

func getScratch(w, h int) *scratch {
	n := (w-1)*h + w*(h-1)
	if n > maxEdges {
		panic("contourmap: grid too large: more than 2^31-1 cell edges")
	}
	s := scratchPool.Get().(*scratch)
	if cap(s.next) < n {
		s.next = make([]int32, n)
		for i := range s.next {
			s.next[i] = noEdge
		}
		s.hasIn = make([]bool, n)
		s.pts = make([]Point, n)
	}
	s.next, s.hasIn, s.pts = s.next[:n], s.hasIn[:n], s.pts[:n]
	return s
}

// putScratch returns s to the pool. Call it only after marchingSquares
// returned normally; after a panic s may be dirty and must be dropped.
func putScratch(s *scratch) {
	if cap(s.next) > maxPooledEdges {
		return
	}
	s.from = s.from[:0]
	scratchPool.Put(s)
}

// invalid reports whether v is NaN or ±Inf. Both are treated as no data.
func invalid(v float64) bool { return v-v != 0 }

func marchingSquares(m *ContourMap, z float64, s *scratch) []Contour {
	w, h := m.W, m.H
	nh := (w - 1) * h

	for bi := range m.blkMin {
		// No cell in a block crosses z unless some sample is above it and
		// some sample is at or below it. NaN samples make neither true.
		if !(m.blkMax[bi] > z && m.blkMin[bi] <= z) {
			continue
		}
		bx, by := bi%m.nbx, bi/m.nbx
		x0, x1 := bx*blockSize, min((bx+1)*blockSize, w-1)
		y0, y1 := by*blockSize, min((by+1)*blockSize, h-1)
		for y := y0; y < y1; y++ {
			// Samples x0..x1 of rows y and y+1; sized so the corner loads
			// below need no bounds checks.
			r0 := m.grid[y*w+x0 : y*w+x1+1]
			r1 := m.grid[(y+1)*w+x0 : (y+1)*w+x1+1]
			r1 = r1[:len(r0)]
			for i := 1; i < len(r0); i++ {
				x := x0 + i - 1
				ul, ur, ll, lr := r0[i-1], r0[i], r1[i-1], r1[i]

				var c int
				if ul > z {
					c |= 1
				}
				if ur > z {
					c |= 2
				}
				if ll > z {
					c |= 4
				}
				if lr > z {
					c |= 8
				}
				if c == 0 || c == 15 {
					continue
				}
				if invalid(ul) || invalid(ur) || invalid(ll) || invalid(lr) {
					continue // no data
				}

				segs := cellSegments[c&15]
				if (c == 6 || c == 9) && connectHigh(m.Saddle, ul, ur, ll, lr, z) {
					segs = saddleHigh[c&15]
				}

				edges := [4]int32{
					eTop:    int32(y*(w-1) + x),
					eBottom: int32((y+1)*(w-1) + x),
					eLeft:   int32(nh + y*w + x),
					eRight:  int32(nh + y*w + x + 1),
				}
				for _, sg := range segs {
					a, b := edges[sg.from&3], edges[sg.to&3]
					s.pts[a] = cellPoint(sg.from, x, y, ul, ur, ll, lr, z)
					s.pts[b] = cellPoint(sg.to, x, y, ul, ur, ll, lr, z)
					s.next[a] = b
					s.hasIn[b] = true
					s.from = append(s.from, a)
				}
			}
		}
	}

	n := len(s.from)
	if n == 0 {
		return nil
	}

	// Open contours start on an edge nothing leaves through: the grid border
	// or the border of a no-data region. Count them to size the output.
	nOpen := 0
	for _, e := range s.from {
		if !s.hasIn[e] {
			nOpen++
		}
	}

	// Each walk emits its start point plus one point per segment. A closed
	// contour has at least 4 segments, so this capacity is rarely exceeded.
	buf := make([]Point, 0, n+nOpen+n/4+1)
	contours := make([]Contour, 0, nOpen+1)

	for _, e := range s.from {
		if !s.hasIn[e] && s.next[e] >= 0 {
			var c Contour
			buf, c = s.walk(e, buf)
			contours = append(contours, c)
		}
	}
	for _, e := range s.from {
		if s.next[e] >= 0 {
			var c Contour
			buf, c = s.walk(e, buf)
			contours = append(contours, c)
		}
	}

	for _, e := range s.from {
		s.next[e] = noEdge
	}
	s.from = s.from[:0]
	return contours
}

// walk follows segments from start until the chain ends or returns to start,
// appending points to buf. It clears the next and hasIn entries it follows.
func (s *scratch) walk(start int32, buf []Point) ([]Point, Contour) {
	i := len(buf)
	buf = append(buf, s.pts[start])
	for e := start; ; {
		n := s.next[e]
		if n < 0 {
			break
		}
		s.next[e] = consumed
		s.hasIn[n] = false
		buf = append(buf, s.pts[n])
		if n == start {
			break
		}
		e = n
	}
	return buf, Contour(buf[i:len(buf):len(buf)])
}
