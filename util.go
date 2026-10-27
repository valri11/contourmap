package contourmap

import (
	"image"
	"image/draw"
)

func imageToGray16(im image.Image) *image.Gray16 {
	b := im.Bounds()
	dst := image.NewGray16(b)
	draw.Draw(dst, b, im, b.Min, draw.Src)
	return dst
}
