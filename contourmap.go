package contourmap

import (
	"image"
	"math"
	"sort"
)

// SaddleRule selects how a saddle cell is resolved. A saddle cell has its two
// high corners on one diagonal and its two low corners on the other, so the
// four samples alone cannot tell whether the high corners connect.
type SaddleRule int

const (
	// SaddleConnectLow never joins the high corners inside a cell. This is
	// the upstream fogleman/contourmap behaviour and the zero value.
	SaddleConnectLow SaddleRule = iota
	// SaddleConnectHigh always joins the high corners inside a cell.
	SaddleConnectHigh
	// SaddleCenter joins the high corners when the mean of the four corners
	// is above the level.
	SaddleCenter
	// SaddleAsymptotic joins the high corners when the saddle point of the
	// bilinear surface through the four corners is above the level. This is
	// exact for a bilinear interpolant.
	SaddleAsymptotic
)

type ContourMap struct {
	W   int     // width of the contour map in pixels
	H   int     // height of the contour map in pixels
	Min float64 // minimum value contained in this contour map
	Max float64 // maximum value contained in this contour map

	// Saddle selects how ambiguous saddle cells are resolved. It must not be
	// changed while Contours is running on another goroutine.
	Saddle SaddleRule

	grid []float64

	// Cells are grouped into blockSize x blockSize blocks. blkMin and blkMax
	// bound the samples each block touches, so a block that cannot cross a
	// level is skipped without classifying its cells.
	nbx    int
	blkMin []float64
	blkMax []float64
}

const blockSize = 16

// FromFloat64s returns a new ContourMap for the provided 2D grid of values.
// len(grid) must equal w * h. NaN and ±Inf values are treated as no data:
// cells that touch one produce no contour segments, so contours end at the
// edge of a no-data region. -math.MaxFloat64 is reserved for the border
// added by Closed and must not appear in the data.
//
// The map keeps a reference to grid and caches per-block value ranges, so
// grid must not be modified after this call.
func FromFloat64s(w, h int, grid []float64) *ContourMap {
	m := &ContourMap{W: w, H: h, Min: math.Inf(1), Max: math.Inf(-1), grid: grid}
	for _, v := range grid {
		if invalid(v) || v == closed {
			continue
		}
		if v < m.Min {
			m.Min = v
		}
		if v > m.Max {
			m.Max = v
		}
	}
	if w < 2 || h < 2 {
		return m
	}

	m.nbx = (w - 2 + blockSize) / blockSize
	nby := (h - 2 + blockSize) / blockSize
	m.blkMin = make([]float64, m.nbx*nby)
	m.blkMax = make([]float64, m.nbx*nby)
	for by := 0; by < nby; by++ {
		y0, y1 := by*blockSize, min((by+1)*blockSize, h-1)
		for bx := 0; bx < m.nbx; bx++ {
			x0, x1 := bx*blockSize, min((bx+1)*blockSize, w-1)
			lo, hi := math.Inf(1), math.Inf(-1)
			for y := y0; y <= y1; y++ {
				for _, v := range grid[y*w+x0 : y*w+x1+1] {
					if v < lo {
						lo = v
					}
					if v > hi {
						hi = v
					}
				}
			}
			m.blkMin[by*m.nbx+bx], m.blkMax[by*m.nbx+bx] = lo, hi
		}
	}
	return m
}

// FromFunction returns a new ContourMap for the provided function.
// The function will be called for all points x = [0, w) and y = [0, h) to
// determine the Z value at each point.
func FromFunction(w, h int, f Function) *ContourMap {
	grid := make([]float64, w*h)
	i := 0
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			grid[i] = f(x, y)
			i++
		}
	}
	return FromFloat64s(w, h, grid)
}

// FromImage returns a new ContourMap for the provided image. The image is
// converted to 16-bit grayscale and will have Z values mapped from
// [0, 65535] to [0, 1].
func FromImage(im image.Image) *ContourMap {
	gray := imageToGray16(im)
	w := gray.Bounds().Size().X
	h := gray.Bounds().Size().Y
	grid := make([]float64, w*h)
	j := 0
	for i := range grid {
		x := int(gray.Pix[j])<<8 | int(gray.Pix[j+1])
		grid[i] = float64(x) / 0xffff
		j += 2
	}
	return FromFloat64s(w, h, grid)
}

func (m *ContourMap) at(x, y int) float64 {
	return m.grid[y*m.W+x]
}

// HistogramZs returns numLevels Z values that split the grid's valid samples
// into equally populated bands. NaN and padding added by Closed are ignored.
func (m *ContourMap) HistogramZs(numLevels int) []float64 {
	result := make([]float64, numLevels)

	values := make([]float64, 0, len(m.grid))
	for _, v := range m.grid {
		if v == v && v != closed {
			values = append(values, v)
		}
	}
	if len(values) == 0 {
		return result
	}
	sort.Float64s(values)

	for i := 0; i < numLevels; i++ {
		t := (float64(i) + 0.5) / float64(numLevels)
		idx := int(t*float64(len(values))) - 1
		if idx < 0 {
			idx = 0
		}
		result[i] = values[idx]
	}
	return result
}

// Contours returns a list of contours that represent isolines at the specified
// Z value. Output is deterministic: the same map and level always produce the
// same contours, in the same order, starting at the same points.
//
// Contours is safe for concurrent use.
func (m *ContourMap) Contours(z float64) []Contour {
	if m.W < 2 || m.H < 2 {
		return nil
	}
	s := getScratch(m.W, m.H)
	out := marchingSquares(m, z, s)
	putScratch(s) // not deferred: a scratch abandoned by a panic is dropped
	return out
}

// ContourLevels returns the contours for each level in zs. It is equivalent
// to calling Contours for each level, but reuses internal buffers across
// levels.
func (m *ContourMap) ContourLevels(zs []float64) [][]Contour {
	out := make([][]Contour, len(zs))
	if m.W < 2 || m.H < 2 {
		return out
	}
	s := getScratch(m.W, m.H)
	for i, z := range zs {
		out[i] = marchingSquares(m, z, s)
	}
	putScratch(s)
	return out
}

// Closed returns a new ContourMap that will ensure all Contours are closed
// paths by following the border when they would normally stop at the edge
// of the grid. Points from the returned map are offset by (1, 1) relative to
// the original grid because of the added border.
//
// The border is below every value, so a level below Min yields one ring
// around the whole grid, which is the base band of a filled contour map.
// Contours that end at a no-data region inside the grid stay open.
func (m *ContourMap) Closed() *ContourMap {
	w := m.W + 2
	h := m.H + 2
	grid := make([]float64, w*h)
	for i := range grid {
		grid[i] = closed
	}
	for y := 0; y < m.H; y++ {
		i := (y+1)*w + 1
		j := y * m.W
		copy(grid[i:], m.grid[j:j+m.W])
	}
	c := FromFloat64s(w, h, grid)
	c.Saddle = m.Saddle
	return c
}
