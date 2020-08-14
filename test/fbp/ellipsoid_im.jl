# ellipsoid_im.jl

using MIRT: ellipsoid_im, shepp_logan_3d_parameters
using MIRT: image_geom, jim
import MIRT: _ellipsoid_im_check_fov

using Plots
using Test: @test, @test_throws


"""
ellipsoid_im_show()
"""
function ellipsoid_im_show(ig)
	x1 = ellipsoid_im(ig, :spheroid)
	p1 = jim(x1, title="spheroid")

	x2 = ellipsoid_im(ig, :zhu; hu_scale=1000)
	p2 = jim(x2, title="zhu", clim=(900,1100))

	plot(p1, p2)
end


ig = image_geom(nx=512, nz=64*2, dz=0.625, fov=500)
ig = ig.down(8)

#test different ellipses
#@inferred # todo
ellipsoid_im(ig)
ellipsoid_im(ig, :zhu)
ellipsoid_im(ig, :kak, oversample=2)
ellipsoid_im(ig, :spheroid, return_params=true)
ellipsoid_im(ig; checkfov=true)
ellipsoid_im(ig; showmem=true)

@test_throws String ellipsoid_im(ig, :bad)
@test_throws String ellipsoid_im(ig, :spheroid, how=:bad)
@test_throws String shepp_logan_3d_parameters(0, 0, 0, :bad)

tmp = [1, 1, 1, 0, 0, 0]
@test !_ellipsoid_im_check_fov(2, 9, 9, [3 0 0 1 1 1 0 0 1], tmp...)
@test !_ellipsoid_im_check_fov(9, 2, 9, [0 3 0 1 1 1 0 0 1], tmp...)
@test !_ellipsoid_im_check_fov(9, 9, 2, [0 0 3 1 1 1 0 0 1], tmp...)

#ell1 = ellipsoid_im(ig, :zhu; how=:fast) # fast doesn't work
#ell2 = ellipsoid_im(ig, :zhu; how=:lowmem) # lowmem calls fast - doesn't work
# todo: test slow vs fast

ellipsoid_im_show(ig)
