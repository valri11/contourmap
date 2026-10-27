## contourmap

Compute contour lines (isolines) for any 2D data in Go.

This is a maintained fork of [fogleman/contourmap](https://github.com/fogleman/contourmap).
It fixes open contours losing their last segment, makes output deterministic,
adds no-data (NaN, ±Inf) support and selectable saddle rules, and is about 7× faster
with about 20× less garbage on tile-sized grids. See
[docs/design/marching-squares.md](docs/design/marching-squares.md) for how it
works and how it differs from upstream, or open
[docs/design/marching-squares.html](docs/design/marching-squares.html) for an
interactive walkthrough.

### Installation

    go get github.com/valri11/contourmap

### Documentation

https://pkg.go.dev/github.com/valri11/contourmap

### Example Usage

#### Creating a ContourMap

A new `ContourMap` can be generated in many different ways, depending on what type of data you have.

Use `FromFloat64s` if you have an array of numbers. The length of the array must equal `width * height`. The two-dimensional data is stored in a flat array in row-major order. NaN values mean no data: contours stop at the edge of a NaN region.

```go
m := contourmap.FromFloat64s(width, height, data)
```

Use `FromImage` if you have an `image.Image`, such as a grayscale heightmap.

```go
m := contourmap.FromImage(im)
```

Use `FromFunction` to specify an arbitrary function that will provide a Z for any given X, Y coordinate.
The function will be called for all points `x = [0, w)` and `y = [0, h)` to determine the Z value at each point in the grid.

```go
var f func(x, y int) float64
...
m := contourmap.FromFunction(width, height, f)
```

#### Finding Contour Lines

Once your `ContourMap` is created, you can use the `Contours` function to find isolines at any given Z height. This function returns a list of contours where each contour is a list of X, Y points.
A `Contour` may be open or closed. Closed contours have `c[0] == c[len(c)-1]`.

```go
contours := m.Contours(z)
for _, contour := range contours {
    for _, point := range contour {
        // do something with points...
        fmt.Println(point.X, point.Y)
    }
}
```

For many levels at once, `ContourLevels` reuses its working buffers across levels:

```go
for i, contours := range m.ContourLevels(levels) {
    // contours for levels[i]
}
```

Output is deterministic, and both methods are safe for concurrent use.

#### Resolving Saddles

A saddle cell has two high corners on one diagonal, so it can be drawn with the high corners joined or kept apart. The default, `SaddleConnectLow`, matches upstream. `SaddleAsymptotic` is exact for the bilinear surface through the cell's corners.

```go
m.Saddle = contourmap.SaddleAsymptotic // or SaddleCenter, SaddleConnectHigh
```

#### Closing Contours at the Grid Perimeter

Contours may end at the edge of the grid data, forming open contours. If you want to force all contours to be closed by following the perimeter of the grid, you can use `ContourMap.Closed` which will generate a new ContourMap that can be used for this purpose. Points from the closed map are offset by (1, 1) because of the added border.

```go
m = m.Closed() // now all contours will be closed paths
```

### Examples

The examples are a separate module, so their drawing dependencies are not pulled into the library.

    $ cd examples
    $ go run ./iceland iceland/iceland.jpg

![Iceland Example](https://i.imgur.com/fd7fUnt.png)

    $ go run ./function

![Function Example](https://i.imgur.com/lbGXPC9.png)

### Testing and benchmarks

    $ go test -race ./
    $ go test -run '^$' -bench . -benchmem ./

The upstream implementation is kept in `legacy_test.go` as a reference for equivalence tests and benchmarks.

### License

MIT. See [LICENSE.md](LICENSE.md).
