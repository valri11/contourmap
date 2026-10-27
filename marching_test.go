package contourmap

import (
	"fmt"
	"math"
	"math/rand"
	"reflect"
	"sort"
	"sync"
	"testing"
)

// terrain returns a deterministic w*h heightfield made of gaussian bumps,
// a ripple and a little noise, similar in character to a DTM tile.
func terrain(w, h int, seed int64) []float64 {
	r := rand.New(rand.NewSource(seed))
	type bump struct{ x, y, s, a float64 }
	bumps := make([]bump, 16)
	for i := range bumps {
		a := 50 + r.Float64()*210
		if r.Float64() < 0.22 {
			a = -a
		}
		bumps[i] = bump{r.Float64() * float64(w), r.Float64() * float64(h), (0.05 + r.Float64()*0.16) * float64(w), a}
	}
	g := make([]float64, w*h)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			v := 140.0
			for _, b := range bumps {
				dx, dy := float64(x)-b.x, float64(y)-b.y
				v += b.a * math.Exp(-(dx*dx+dy*dy)/(2*b.s*b.s))
			}
			v += 14*math.Sin(float64(x)*0.11+2*math.Cos(float64(y)*0.05))*math.Cos(float64(y)*0.08) + r.Float64()*0.5
			g[y*w+x] = v
		}
	}
	return g
}

func levels(m *ContourMap, interval float64) []float64 {
	var zs []float64
	for k := math.Ceil(m.Min / interval); k*interval <= m.Max; k++ {
		zs = append(zs, k*interval-1e-4)
	}
	return zs
}

func TestOpenContourKeepsEndpoints(t *testing.T) {
	m := FromFloat64s(3, 3, []float64{0, 1, 0, 0, 0, 0, 0, 0, 0})
	got := m.Contours(0.5)
	want := []Contour{{{1.5, 0}, {1, 0.5}, {0.5, 0}}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("got %v, want %v", got, want)
	}

	// upstream drops the final point of every open contour
	legacy := legacyMarchingSquares(m, m.W, m.H, 0.5)
	if len(legacy) != 1 || len(legacy[0]) != 2 {
		t.Fatalf("legacy behaviour changed: %v", legacy)
	}
}

func TestPointsLieOnLevel(t *testing.T) {
	const w, h = 64, 48
	g := terrain(w, h, 1)
	m := FromFloat64s(w, h, g)
	for _, z := range levels(m, 10) {
		for _, c := range m.Contours(z) {
			for _, p := range c {
				var a, b, f float64
				switch {
				case p.X == math.Trunc(p.X):
					x, y := int(p.X), int(p.Y)
					a, b, f = m.at(x, y), m.at(x, y+1), p.Y-float64(y)
				case p.Y == math.Trunc(p.Y):
					x, y := int(p.X), int(p.Y)
					a, b, f = m.at(x, y), m.at(x+1, y), p.X-float64(x)
				default:
					t.Fatalf("point %v is not on a grid edge", p)
				}
				if v := a + (b-a)*f; math.Abs(v-z) > 1e-6 {
					t.Fatalf("point %v interpolates to %v, want %v", p, v, z)
				}
			}
		}
	}
}

func segmentKeys(cs []Contour) map[string]int {
	keys := map[string]int{}
	for _, c := range cs {
		for i := 1; i < len(c); i++ {
			keys[fmt.Sprintf("%.9f,%.9f>%.9f,%.9f", c[i-1].X, c[i-1].Y, c[i].X, c[i].Y)]++
		}
	}
	return keys
}

func isClosed(c Contour) bool { return len(c) > 2 && c[0] == c[len(c)-1] }

func TestMatchesLegacyPlusMissingEndSegments(t *testing.T) {
	for seed := int64(1); seed <= 5; seed++ {
		const w, h = 80, 60
		m := FromFloat64s(w, h, terrain(w, h, seed))
		for _, z := range levels(m, 15) {
			got := m.Contours(z)
			old := legacyMarchingSquares(m, w, h, z)

			gk, ok := segmentKeys(got), segmentKeys(old)
			for k, n := range ok {
				if gk[k] != n {
					t.Fatalf("seed %d z %v: legacy segment %s missing", seed, z, k)
				}
			}
			open := 0
			for _, c := range got {
				if !isClosed(c) {
					open++
				}
			}
			if len(gk) != len(ok)+open || len(got) != len(old) {
				t.Fatalf("seed %d z %v: %d segments / %d contours, legacy %d / %d, open %d",
					seed, z, len(gk), len(got), len(ok), len(old), open)
			}
		}
	}
}

func TestDeterministicAndConcurrent(t *testing.T) {
	const w, h = 100, 70
	m := FromFloat64s(w, h, terrain(w, h, 3))
	zs := levels(m, 8)
	want := m.ContourLevels(zs)

	var wg sync.WaitGroup
	for i := 0; i < 8; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j, z := range zs {
				if got := m.Contours(z); !reflect.DeepEqual(got, want[j]) {
					t.Errorf("level %v: output differs between calls", z)
					return
				}
			}
		}()
	}
	wg.Wait()
}

func TestClosedProducesOnlyClosedContours(t *testing.T) {
	const w, h = 50, 40
	m := FromFloat64s(w, h, terrain(w, h, 4)).Closed()
	for _, z := range levels(m, 12) {
		for _, c := range m.Contours(z) {
			if !isClosed(c) {
				t.Fatalf("z %v: open contour %v", z, c)
			}
		}
	}
}

func TestSaddleRules(t *testing.T) {
	// Two peaks of 9 on a diagonal; the centre cell is case 9 with its
	// bilinear saddle at 5.5.
	g := []float64{
		0, 0, 0, 0, 0, 0,
		0, 3, 3, 1, 0, 0,
		0, 3, 9, 2, 1, 0,
		0, 1, 2, 9, 3, 0,
		0, 0, 1, 3, 3, 0,
		0, 0, 0, 0, 0, 0,
	}
	cases := []struct {
		rule SaddleRule
		z    float64
		want int
	}{
		{SaddleConnectLow, 5.2, 2},
		{SaddleConnectHigh, 5.8, 1},
		{SaddleCenter, 5.2, 1},
		{SaddleCenter, 5.8, 2},
		{SaddleAsymptotic, 5.2, 1},
		{SaddleAsymptotic, 5.8, 2},
	}
	for _, tc := range cases {
		m := FromFloat64s(6, 6, g)
		m.Saddle = tc.rule
		if got := len(m.Contours(tc.z)); got != tc.want {
			t.Errorf("rule %d z %v: %d contours, want %d", tc.rule, tc.z, got, tc.want)
		}
	}
}

func TestNaNIsNoData(t *testing.T) {
	const w, h = 40, 30
	g := terrain(w, h, 5)
	for y := 10; y < 18; y++ {
		for x := 12; x < 25; x++ {
			g[y*w+x] = math.NaN()
		}
	}
	m := FromFloat64s(w, h, g)
	if math.IsNaN(m.Min) || math.IsNaN(m.Max) {
		t.Fatalf("Min/Max poisoned by NaN: %v %v", m.Min, m.Max)
	}
	for _, z := range levels(m, 10) {
		for _, c := range m.Contours(z) {
			for _, p := range c {
				if math.IsNaN(p.X) || math.IsNaN(p.Y) {
					t.Fatalf("NaN point in %v", c)
				}
				if p.X > 12 && p.X < 24 && p.Y > 10 && p.Y < 17 {
					t.Fatalf("point %v inside no-data hole", p)
				}
			}
		}
	}
}

func TestHistogramZsMatchesUpstream(t *testing.T) {
	const w, h = 30, 20
	m := FromFloat64s(w, h, terrain(w, h, 6))

	// upstream algorithm
	hist := map[float64]int{}
	for _, v := range m.grid {
		hist[v]++
	}
	keys := make([]float64, 0, len(hist))
	for k := range hist {
		keys = append(keys, k)
	}
	sort.Float64s(keys)
	const n = 12
	want := make([]float64, n)
	for i := 0; i < n; i++ {
		count := int((float64(i) + 0.5) / n * float64(len(m.grid)))
		total := 0
		for _, k := range keys {
			total += hist[k]
			if total >= count {
				want[i] = k
				break
			}
		}
	}
	if got := m.HistogramZs(n); !reflect.DeepEqual(got, want) {
		t.Fatalf("got %v, want %v", got, want)
	}
}

func TestScratchIsCleanAfterUse(t *testing.T) {
	const w, h = 60, 40
	m := FromFloat64s(w, h, terrain(w, h, 7))
	s := getScratch(w, h)
	for _, z := range levels(m, 5) {
		marchingSquares(m, z, s)
	}
	for i := range s.next {
		if s.next[i] != noEdge || s.hasIn[i] {
			t.Fatalf("edge %d left dirty: next=%d hasIn=%v", i, s.next[i], s.hasIn[i])
		}
	}
	putScratch(s)
}

func TestDegenerateSizes(t *testing.T) {
	for _, sz := range [][2]int{{0, 0}, {1, 1}, {1, 5}, {5, 1}} {
		m := FromFloat64s(sz[0], sz[1], make([]float64, sz[0]*sz[1]))
		if c := m.Contours(0); c != nil {
			t.Errorf("%v: got %v", sz, c)
		}
	}
}

// The elevation service's shape: a 262x262 window (256 px tile plus 3 px halo)
// contoured at every 10 m between min and max.
func benchSetup() (*ContourMap, []float64) {
	const w, h = 262, 262
	m := FromFloat64s(w, h, terrain(w, h, 42))
	return m, levels(m, 10)
}

func BenchmarkLegacy(b *testing.B) {
	m, zs := benchSetup()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		for _, z := range zs {
			legacyMarchingSquares(m, m.W, m.H, z)
		}
	}
}

func BenchmarkContours(b *testing.B) {
	m, zs := benchSetup()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		for _, z := range zs {
			m.Contours(z)
		}
	}
}

func BenchmarkContourLevels(b *testing.B) {
	m, zs := benchSetup()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		m.ContourLevels(zs)
	}
}

// New map per request, as the service builds one per tile.
func BenchmarkPerRequest(b *testing.B) {
	m, zs := benchSetup()
	grid := m.grid
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		FromFloat64s(262, 262, grid).ContourLevels(zs)
	}
}

func TestInfIsNoData(t *testing.T) {
	g := terrain(20, 20, 8)
	g[5*20+5], g[12*20+9] = math.Inf(1), math.Inf(-1)
	m := FromFloat64s(20, 20, g)
	if math.IsInf(m.Min, 0) || math.IsInf(m.Max, 0) {
		t.Fatalf("Min/Max include Inf: %v %v", m.Min, m.Max)
	}
	for _, z := range levels(m, 5) {
		for _, c := range m.Contours(z) {
			for _, p := range c {
				if math.IsNaN(p.X) || math.IsNaN(p.Y) {
					t.Fatalf("NaN point in %v", c)
				}
			}
		}
	}
}

func TestLargeScratchIsNotPooled(t *testing.T) {
	s := &scratch{next: make([]int32, 0, maxPooledEdges+1)}
	putScratch(s)
	for i := 0; i < 4; i++ {
		if got := scratchPool.Get().(*scratch); got == s {
			t.Fatal("oversized scratch was pooled")
		}
	}
}
