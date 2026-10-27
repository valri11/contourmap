# Marching squares in contourmap

This document explains how contourmap turns a 2D grid of values into contour
lines (isolines), and why the implementation is shaped the way it is.

For an interactive version with live demos of every step, open
[marching-squares.html](marching-squares.html) in a browser. It is a single
self-contained file; the only network request is for web fonts.

## Contents

1. [The problem](#1-the-problem)
2. [Cells and corner classification](#2-cells-and-corner-classification)
3. [The case table](#3-the-case-table)
4. [Edge interpolation](#4-edge-interpolation)
5. [Saddles](#5-saddles)
6. [Stitching segments into polylines](#6-stitching-segments-into-polylines)
7. [Levels, no-data and closed contours](#7-levels-no-data-and-closed-contours)
8. [Tiled grids](#8-tiled-grids)
9. [Performance design](#9-performance-design)
10. [Differences from upstream](#10-differences-from-upstream)

## 1. The problem

The input is a `w × h` grid of samples in row-major order: `grid[y*w + x]`.
A contour at level `z` is the curve where the underlying surface crosses `z`.
Only the samples are known, so the algorithm has to estimate where the curve
passes between them.

Coordinates in the output use grid units: sample `(x, y)` is at point
`(x, y)`, and a crossing between two samples lies on the edge joining them.

## 2. Cells and corner classification

Marching squares treats each 2 × 2 block of neighbouring samples as a *cell*.
A `w × h` grid has `(w−1) × (h−1)` cells, and each one is handled on its own.

Each corner gets one bit: 1 if its value is strictly greater than `z`. The
bits are packed as:

```
ul = 1   ur = 2
ll = 4   lr = 8
```

which gives a case number from 0 to 15. Cases 0 (all low) and 15 (all high)
contain no crossing and are skipped. On real terrain most cells are one of
these two at any given level.

The comparison is strict (`v > z`). When samples hold exact integers and the
level is also an integer, a sample equal to the level sits exactly on the line.
That produces zero-length segments, so callers usually contour at `z − ε`
(for example `z − 1e-4`).

## 3. The case table

An edge of the cell is crossed exactly when its two endpoints disagree. The
case number fixes which edges are crossed and how to pair them into segments.
Edges are named by position: `t` top, `b` bottom, `l` left, `r` right.

| Case | Bits (lr ll ur ul) | High corners | Segments (entry → exit) |
|-----:|:------------------:|--------------|-------------------------|
| 0    | 0000 | none        | none |
| 1    | 0001 | ul          | t → l |
| 2    | 0010 | ur          | r → t |
| 3    | 0011 | ul ur       | r → l |
| 4    | 0100 | ll          | l → b |
| 5    | 0101 | ul ll       | t → b |
| 6    | 0110 | ur ll       | r → t, l → b *(saddle)* |
| 7    | 0111 | ul ur ll    | r → b |
| 8    | 1000 | lr          | b → r |
| 9    | 1001 | ul lr       | t → l, b → r *(saddle)* |
| 10   | 1010 | ur lr       | b → t |
| 11   | 1011 | ul ur lr    | b → l |
| 12   | 1100 | ll lr       | l → r |
| 13   | 1101 | ul ll lr    | t → r |
| 14   | 1110 | ur ll lr    | l → t |
| 15   | 1111 | all         | none |

Segments are **directed**: walking from entry to exit always keeps the high
side on the same hand. This is what makes stitching cheap (§6). The table
lives in `cellSegments` in `marching.go`.

## 4. Edge interpolation

On a crossed edge with endpoint values `a` and `b`, the crossing is placed by
linear interpolation:

```
t = (z − a) / (b − a)        // 0 at a, 1 at b
```

`t` is clamped to `[1e-9, 1 − 1e-9]` so a crossing never lands exactly on a
sample. Without interpolation every crossing would sit at the midpoint of its
edge and the lines would come out jagged.

Two neighbouring cells compute the point on their shared edge from the same
two samples in the same order, so they agree on it bit for bit.

## 5. Saddles

Cases 6 and 9 have the two high corners on one diagonal and the two low
corners on the other. Both pairings in the table and the alternative are
valid, and the four samples cannot say which one is right: either the high
corners connect through the middle of the cell (one ridge) or they don't (two
hills).

`ContourMap.Saddle` selects the rule:

| Rule | High corners joined when | Notes |
|------|--------------------------|-------|
| `SaddleConnectLow` (default) | never | Upstream behaviour. Consistent, sometimes topologically wrong. |
| `SaddleConnectHigh` | always | The mirror image of the default. |
| `SaddleCenter` | `(ul + ur + ll + lr) / 4 > z` | Cheap approximation. |
| `SaddleAsymptotic` | `(ul·lr − ur·ll) / (ul + lr − ur − ll) > z` | Saddle value of the bilinear surface through the four corners. Exact for that surface. |

The alternative pairings are in `saddleHigh`:
case 6 becomes `l → t, r → b` and case 9 becomes `t → r, b → l`.

Saddle rules matter most for small intervals on flat, noisy data, where
saddles cluster and the fixed rule can split one contour into many small loops.

```go
m := contourmap.FromFloat64s(w, h, grid)
m.Saddle = contourmap.SaddleAsymptotic
```

## 6. Stitching segments into polylines

Per-cell segments have to be joined into long polylines for rendering,
simplification and labelling.

### Edge ids

Every grid edge has an integer id:

```
horizontal edge (x, y)–(x+1, y):  y·(w−1) + x                 x ∈ [0, w−1), y ∈ [0, h)
vertical   edge (x, y)–(x, y+1):  (w−1)·h + y·w + x           x ∈ [0, w),   y ∈ [0, h−1)
```

For the cell whose upper-left sample is `(x, y)`, its four edges are
`t = h(x, y)`, `b = h(x, y+1)`, `l = v(x, y)` and `r = v(x+1, y)`. Neighbouring
cells derive the same id for a shared edge, so no float comparison is needed
to match segments.

### The graph

Because segments are directed, every crossed edge is the **entry** of at most
one segment and the **exit** of at most one segment. The scan fills three
arrays indexed by edge id:

```go
next[a]  = b     // the segment entering through a leaves through b
hasIn[b] = true  // some segment leaves through b
pts[a], pts[b]   // crossing points on both edges
from = append(from, a)   // entry edges, in scan order
```

The result is a set of disjoint chains and cycles.

### Walking

```
for each a in from:                   // open contours
    if !hasIn[a] && next[a] is set:   //   a chain starts where nothing leads in
        walk(a)
for each a in from:                   // whatever remains is a loop
    if next[a] is set:
        walk(a)

walk(s):
    emit pts[s]
    e := s
    loop:
        n := next[e]; if n unset: stop      // open contour ends here
        mark e consumed; hasIn[n] = false
        emit pts[n]
        if n == s: stop                     // closed: first point == last point
        e = n
```

An open contour starts on an edge that nothing exits through: the grid border,
or the border of a no-data region. Finding starts by in-degree handles both
without special cases. A closed contour is emitted with `c[0] == c[len(c)-1]`.

Starts are taken in scan order, so output is deterministic: the same grid and
level always give the same contours, in the same order, from the same start
points.

## 7. Levels, no-data and closed contours

**Levels.** A contour map runs the single-level algorithm once per level.
`ContourLevels(zs)` does this with one set of working buffers. Generate levels
from an integer index, `z = k·interval`, rather than adding `interval` to a
float in a loop, which drifts.

**No-data.** NaN and ±Inf mean no data. `Min` and `Max` ignore them, and a cell with
any such corner emits nothing, so contours stop at the edge of the hole as
open contours. Convert sentinels such as −9999 to NaN before building the
map. Replacing them with 0 creates an artificial cliff, and every level
between the terrain and 0 then runs along it.

**Closed contours.** `Closed()` returns a copy of the map padded with a
one-sample border holding an internal sentinel (`−MaxFloat64`). The border is
always low, so contours that reach it close by following the grid's edge, which is
what filled contour bands need. A level below `Min` gives one ring around the whole grid (the base band), and contours that end at an interior no-data region stay open. Output points are offset by `(1, 1)` from the
original grid. Saddle rules other than the default treat any cell touching
the border as "not connected".

## 8. Tiled grids

Map clients request contours one tile at a time. If each tile is contoured from
its own samples only, the cells straddling tile borders are never processed:
the last sample of one tile and the first of the next are neighbours, but no
single tile contains both. Every line breaks at every seam.

The fix is a **halo**: contour a window that includes a few samples from each
neighbouring tile, then clip back to the tile. One sample is enough to close
the seams; a few more help when the output is later simplified.

```go
// window: tile plus `halo` samples from each neighbour, row-major
m := contourmap.FromFloat64s(tile+2*halo, tile+2*halo, window)
for _, lines := range m.ContourLevels(levels) {
    for _, c := range lines {
        for _, p := range c {
            px, py := p.X-halo, p.Y-halo // tile pixel coordinates
            _, _ = px, py
        }
    }
}
```

Clip to the tile plus any vector-tile buffer afterwards. Lines inside the
halo duplicate the neighbour tile's lines.

## 9. Performance design

The hot path is classifying cells. Everything else is arranged so that the
classification loop is all the work there is and nothing allocates per cell.

**Block skipping.** `FromFloat64s` groups cells into 16 × 16 blocks and records
the min and max of the samples each block touches. For each level, a block is
classified only if `blockMax > z && blockMin ≤ z`. Most blocks fail that test
at any given level, so the cost per level scales with how much of the grid the
level actually crosses.

**Dense, pooled scratch.** Stitching uses flat slices indexed by edge id
instead of maps. They come from a package-level `sync.Pool` and are left in a
clean state at the end of every call (every `next` reset, every `hasIn` cleared
by the walk), so they are reused without being zeroed. A map built per request
for a same-sized grid picks up the warm buffers. Scratch size is about
`2·w·h × 21` bytes: roughly 2.9 MB for a 262 × 262 grid. Scratch for more than 2²¹ edges (a grid of about 1000 × 1000) is not returned to the pool, so one large grid can't pin tens of megabytes for later small calls. A scratch is only returned after a normal return, never after a panic. The map keeps a reference to the input slice and caches block ranges, so the input must not be modified afterwards.

**Output allocation.** All points for one level go into one backing slice whose
capacity is computed before the walk (`segments + openContours + segments/4`,
since a loop has at least four segments). Each `Contour` is a full slice
expression into it (`buf[i:j:j]`), so appending to one contour can never
overwrite the next. That gives roughly two allocations per level.

**Concurrency.** A `ContourMap` is read-only after construction. `Contours` and
`ContourLevels` are safe to call from multiple goroutines; each call takes its
own scratch from the pool. Don't change `Saddle` while calls are running.

**Benchmark** (Apple M2 Pro, Go 1.27, 262 × 262 synthetic terrain, all 41
levels at 10 m intervals; `go test -bench . -benchmem`):

| Implementation | Time | Bytes | Allocations |
|----------------|-----:|------:|------------:|
| upstream (`legacy_test.go`) | 9.7 ms | 12.8 MB | 2,459 |
| this library, `Contours` per level | 1.36 ms | 0.58 MB | 97 |
| this library, new map + `ContourLevels` | 1.57 ms | 0.58 MB | 101 |

## 10. Differences from upstream

This module is derived from
[github.com/fogleman/contourmap](https://github.com/fogleman/contourmap).

| Area | Upstream | This library |
|------|----------|--------------|
| Open contours | Last point dropped: crossing points were stored only for entry edges, so the final segment at the border was lost | Endpoints kept |
| Output order | Random (Go map iteration) | Deterministic (scan order) |
| Stitching storage | Two maps per level with struct keys | Pooled flat slices by edge id |
| Cells scanned | Every cell, every level | Only 16 × 16 blocks whose range spans the level |
| Saddles | Fixed `connectHigh = false` | `Saddle` rule, same default |
| NaN, ±Inf | NaN counted as low; NaN points; `Min`/`Max` poisoned | No data; ignored by `Min`/`Max` |
| `HistogramZs` | Map histogram, O(levels × unique values) | Sort once; ignores NaN and `Closed()` padding |
| `FromImage` | Wrong source point for images whose bounds don't start at (0, 0) | Uses `Bounds().Min` |
| Multi-level API | none | `ContourLevels` |

The upstream implementation is kept in `legacy_test.go`. The tests check that
on random terrain the new output contains every upstream segment, plus exactly
one extra segment per open contour.
