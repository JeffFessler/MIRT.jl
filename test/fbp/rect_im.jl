# rect_im.jl

using MIRT: rect_im, image_geom, jim

using Plots: plot
using Test: @test, @test_throws, @inferred


"""
rect_im_show()
"""
function rect_im_show()
	ig = image_geom(nx=2^8, ny=2^8, fov=100)

	x1 = rect_im(ig, [[[0.5, 0, 3, 20]*ig.dx..., 0, 1]';], oversample=3)
	p1 = jim(x1)

	x2 = rect_im(ig, :default ; oversample=3, chat=true)
	p2 = jim(x2, title="default rects")

	x3 = rect_im(ig, :my_rect ; oversample=3, chat=true)
	p3 = jim(x3, title="my rect")

	x4 =rect_im(ig, :smiley ; oversample=3, chat=true)
	p4 = jim(x4, title="smiley")

	plot(p1, p2, p3, p4)
end



ig = image_geom(nx=2^8, dx=3)
@test_throws String rect_im(ig, :bad)
rect_im(ig, :smiley, how=:fast, replace=true)
rect_im(ig, how=:slow, replace=true)
rect_im(32, ny=30, dx=3, params=:default, how=:slow, return_params=true)
rect_im(32, 30)
@test_throws String rect_im(32, :default, how=:bad)

rect_im_show()
