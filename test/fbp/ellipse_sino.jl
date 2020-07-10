# ellipse_sino.jl

using MIRT: ellipse_sino
using MIRT: image_geom, sino_geom, ellipse_im, jim, prompt
using Plots: plot, Plot
using Test: @test, @inferred


"""
`ellipse_sino_show()`
show examples

To see more, use `MIRT.ellipse_sino_show(ip=2)`
"""
function ellipse_sino_show( ;
	down::Int = 4,
	ell::Matrix = [ 40 60 50 150 20 10 ],
	orbit::Real = 360,
	na::Int = 400,
	oversample::Int = 2,
	ip::Int = 1,
)

	ig = image_geom(nx=512, ny=504, dx=1, dy=1, mask=:circ)
	ig = ig.down(down)
	xtrue = ellipse_im(ig, ell ; oversample=2)
	ig = image_geom(nx=ig.nx, ny=ig.ny, dx=ig.dx, dy=ig.dy, mask = xtrue .> 0)

	geoms = (
		sino_geom(:par, nb = 444, na=na, down=down, d = 1, orbit=orbit,
			offset = 0.25),
		sino_geom(:fan, nb = 888, na=na, d = 1.0, orbit=orbit,
			offset = 0.75, dsd = 949, dod = 408, down=down),
		sino_geom(:fan, nb = 888, na=na, d = 1.0, orbit=orbit,
			offset = 0.75, dsd = 949, dod = 408, down=down,
			dfs = Inf, source_offset = 0.7), # flat fan
		sino_geom(:moj, nb = 630, na=na, down=down, d = 1.0, orbit=orbit,
			offset = 0.25),
	)

	ngeom = length(geoms)
	pl = Array{Plot}(undef, ngeom, 2)

	for ii=1:ngeom
		sg = geoms[ii]
		sino = ellipse_sino(sg, ell; oversample=oversample)
		dfs = sg.how === :fan ? " dfs=$(sg.dfs)" : ""
		pl[ii,1] = jim(sino, title="$(sg.how)$dfs")
		pl[ii,2] = sg.plot(ig=ig)
	end
#	plot(pl..., layout=(2,ngeom)) # too small
	plot(pl[:,ip]...)
end


ell = [ 40 70 50 150 20 10 ]
sg = sino_geom(:ge1, down=8)
#@inferred # todo
ellipse_sino(sg, ell; xscale=-1, yscale=-1) # test scale

ellipse_sino_show()
prompt()
ellipse_sino_show(ip=2)
