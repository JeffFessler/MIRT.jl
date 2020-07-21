# image_geom.jl

using MIRT: ImageGeom, image_geom, cbct
using MIRT: image_geom_ellipse
#using MIRT: jim

using Test: @test, @testset, @test_throws, @inferred


function image_geom_test2(ig::ImageGeom)
	# test 2D functions provided by the constructor
	ig.dim
	ig.dims
	ig.deltas
	ig.offsets
	ig.x
	ig.y
	ig.wx
	ig.wy
	ig.xg
	ig.yg
	ig.fovs
	ig.np
	ig.mask_outline
	ig.ones
	ig.zeros
	ig.u
	ig.v
	ig.ug
	ig.vg
	ig.fg
	@test ig.shape(vec(ig.ones)) == ig.ones
	@test ig.embed(ig.ones[ig.mask]) == Float32.(ig.mask)
	@test ig.maskit(ig.ones) == Float32.(ones(ig.np))
#	@inferred # todo
	ig.unitv()
	ig.unitv(j=4)
	ig.unitv(i=ntuple(i->1, length(ig.dim)))
	ig.unitv(c=ntuple(i->0, length(ig.dim)))
#= todo-i: why do these fail?
	@inferred image_geom_ellipse(8, 10, 1, 2)
	@inferred ig.circ()
	@inferred ig.plot()
	@inferred ig.unitv()
=#
	ig.circ()
	ig.plot()
#	isinteractive() && gui()
	ig.unitv()
	@inferred ig.down(2)
	@inferred ig.over(2)
	image_geom_ellipse(8, 10, 1, 2)
	true
end


function image_geom_test2()
	#@inferred
	image_geom(nx=16, dx=2, offsets=:dsp, mask=:all_but_edge_xy)
	@test_throws String image_geom(nx=16, dx=1, offsets=:bad)
	@test_throws String image_geom(nx=16, dx=1, mask=:bad) # mask type
	@test_throws String image_geom(nx=16, dx=1, mask=trues(2,2)) # mask size
	ig = image_geom(nx=16, dx=2)
	show(isinteractive() ? stdout : devnull, ig)
	show(isinteractive() ? stdout : devnull, MIME("text/plain"), ig)
	image_geom_test2(ig)
	ig = image_geom(nx=16, dx=2, mask=:circ)
	ig.over(2)
	ig.down(3) # test both even and non-even factors
	ig.help
	true
end


function image_geom_test3(ig::ImageGeom)
	@test image_geom_test2(ig)
	ig.wz
	ig.zg
	ig.mask_or
	@inferred ig.expand_nz(2)
	@inferred cbct(ig)
	@test_throws String image_geom(nx=1, dx=1, nz=2, offset_z=-1, offsets=:dsp)
	true
end


function image_geom_test3()
	#@inferred
	image_geom(nx=16, nz=4, dx=2, zfov=1)
	ig = image_geom(nx=16, nz=4, dx=2, dz=3)
	image_geom_test3(ig)
	true
end


@testset "2d" begin
	@test all(ImageGeom{1}((2,), (3,), (0,)).mask)
	@test image_geom_test2()
end
@testset "3d" begin
	@test image_geom_test3()
end
