# cuboid_im.jl

using MIRT: cuboid_im, image_geom, jim

using Plots
using Test: @test, @test_throws, @inferred
#using MIRT: jim, MIRT_image_geom, downsample3, rotate2d


"""
cuboid_im_show()
"""
function cuboid_im_show()
	ig = image_geom(nx=2^6, ny=2^6, nz=2^3, dz=1, fov=240)

	x1 = cuboid_im(ig, :default, how=:exact)
	p1 = jim(x1, title="exact")

	x2 = cuboid_im(ig, :default, how=:sample)
	p2 = jim(x2, title="sample")

	x3 = cuboid_im(ig, :default, how=:lowmem1)
	p3 = jim(x3, title="lowmem")

	x4 = cuboid_im(ig, :rotate)
	p4 = jim(x4, title="rotate") # runs sample

	plot(p1, p2, p3, p4)
end


ig = image_geom(nx=2^3, ny=2^3, nz=2^3, dz=1, fov=240)

@test_throws String cuboid_im(ig, :bad)
@test_throws String cuboid_im(ig, :rotate, how=:exact)
@test_throws String cuboid_im(ig, :rotate, how=:bad)
@test_throws String cuboid_im(ig, [0 0 0 0 0 0 0 1 1])
#@test @inferred cuboid_im(ig, :rotate, return_params=true) # todo: fails
@test cuboid_im(ig, :rotate, return_params=true) isa Tuple

diam = abs.([2*ig.dx 2.7*ig.dy 3.2*ig.dz])
params = [1.4 -0.5 1 diam 0 0 1]
phantom_exact = cuboid_im(ig, params, how=:exact)
dxyz = abs(ig.dx * ig.dy * ig.dz)
vol_true = prod(diam)
vol_phantom = sum(phantom_exact) * dxyz
@test vol_true ≈ vol_phantom

phantom_sample = cuboid_im(ig, :default, how=:sample, oversample=3)
vol_phantom2 = sum(phantom_sample) * dxyz
#@test vol_phantom2 ≈ vol_true # todo: fails (check matlab)

xrl = cuboid_im(ig, :rotate, how=:lowmem1)
xrs = cuboid_im(ig, :rotate, how=:sample)
@test xrl == xrs

#x4 = cuboid_im(ig, :rotate, how=:exact) exact does not support rotation
x3 = cuboid_im(ig, :default, how=:lowmem1)
x5 = cuboid_im(ig)

cuboid_im_show()
