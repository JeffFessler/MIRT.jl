# rect_sino.jl

using MIRT: rect_sino, sino_geom, jim
import MIRT: trapezoid

using Plots: plot
using Test: @test, @test_throws, @inferred


"""
`rect_sino_show()`
show examples
"""
function rect_sino_show( ;
	down::Int = 4,
	rect::Matrix = [70 20 70 30 0 1],
	orbit::Real = 360,
	na::Int = 400,
	oversample::Int = 2,
)

#=
	ig = image_geom(nx=512, ny=504, fov=500)
	ig = ig.down(down)
	xtrue = rect_im(ig, rect; oversample=4)
=#

	arg = (na=na, down=down, orbit=orbit, offset=0.25)
	geoms = (
		sino_geom(:par, nb = 444, d = 1 ; arg...),
		sino_geom(:fan, nb = 888, d = 1 ; arg..., dsd = 949, dod = 408),
		sino_geom(:fan, nb = 888, d = 1 ; arg...,
			dsd = 949, dod = 408, dfs = Inf), # flat fan
		sino_geom(:moj, nb = 666, d = 1 ; arg...),
	)

	ngeom = length(geoms)
	pl = Array{Any}(undef, ngeom)

	for ii=1:ngeom
		sg = geoms[ii]
		over = (sg isa SinoMoj) ? 1 : oversample
		sino = rect_sino(sg, rect; oversample=over)
		pl[ii] = jim(sg.r, sg.ad, sino, title="$(typeof(sg))")
	end
	plot(pl...)
end


param = [60 10.5 38 27 0 1]
sg = sino_geom(:ge1, down=8)
#@inferred # todo
rect_sino(sg, param; xscale=-1, yscale=-1) # test scale
@inferred trapezoid(2.5, 1, 2, 3, 4.)
t = LinRange(-2, 5, 101)
f = trapezoid.(t, 1, 2, 3, 4.)

rect_sino_show()
