# sino_geom.jl

using MIRT: sino_geom
using MIRT: sino_geom_help, sino_geom_plot_grids
using MIRT: image_geom, prompt
import MIRT: sino_geom_gamma_dfs

using Plots: plot!, plot, scatter, gui
using Test: @test, @test_throws, @inferred


"""
    sino_geom_show()
show an example of each of the 4 main geometries
"""
function sino_geom_show( ; kwarg...)
	down = 4
	ig = image_geom(nx=512, fov=500)
	ig = ig.down(down)

	arg = (nb = 888, na = down*8)
	nb_moj = round(Int, ig.nx*down*sqrt(2)) # match rfov
	sg_list = (
		sino_geom(:par ; d=ig.fovs[1]/arg.nb, arg..., kwarg...),
		sino_geom(:ge1 ; dfs=0, arg..., kwarg...), # arc
		sino_geom(:ge1 ; dfs=Inf, arg..., nb=960, kwarg...), # flat
		sino_geom(:moj ; d=1, arg..., nb=nb_moj, kwarg...),
		)

	nsg = length(sg_list)
	pl = Array{Any}(undef, nsg)
	for ii = 1:nsg
		sg = sg_list[ii].down(down)
		pl[ii] = plot(); sino_geom_plot!(sg, plot!)
		plot!(pl[ii], xlim=[-1,1]*550, ylim=[-1,1]*550)
		plot!(pl[ii], xtick=[-1,0,1]*250)
	end
	return pl
end


down = 8
ig = image_geom(nx=512, fov=500).down(down)
ig = image_geom(nx=ig.nx, dx=ig.dx, mask=ig.circ())

sg_list = (
	sino_geom(:par ; down=down, d=4),
	sino_geom(:ge1 ; down=down, orbit_start=20, dfs=0), # arc
	sino_geom(:ge1 ; down=down, orbit_start=20, dfs=Inf), # flat
	sino_geom(:moj ; down=down, d=4*sqrt(2)),
	sino_geom(:fan ; down=down, orbit=:short),
)

sg_list[1].help

nsg = length(sg_list)
pl = Array{Any}(undef, nsg)
for ii = 1:nsg
	sg = sg_list[ii]

	sg.ad[2]
	sg.rfov
	#@inferred
	sg.down(2)
	#@inferred
	sg.over(2)
	sg.dim
	sg.w
	sg.ones
	sg.zeros
	sg.dr
	sg.ds
	sg.r
	sg.s
	sg.ad
	sg.ar
	sg.xds
	sg.yds

	if sg isa SinoFan
		sg.gamma
		sg.gamma_max
		sg.dfs
		sg.dso
		sg.orbit_short
	end

	sg.grid
	sg.plot_grid(plot)
#	gui(); prompt()

	if sg isa SinoMoj
		sg.d_moj(0)
		sg.d_ang # angular dependent d for :moj
	end

	@test sg.shape(vec(sg.ones)) == sg.ones
	sg.taufun(ig.x, 0*ig.x)
	sg.unitv()
	sg.unitv(ib=1, ia=2)

	pl[ii] = plot(); sg.plot!(plot! ; ig=ig)
#	gui(); prompt()
end

#@inferred todo
@test sino_geom(:ge1 ; units=:cm) isa SinoFanArc

@test_throws String sino_geom(:badhow)
@test_throws String sino_geom(:ge1 ; dfs=-1)
@test_throws String sino_geom(:ge1 ; units=:bad)

sg = sino_geom(:ge1 ; orbit=:short)
@test sino_geom_gamma_dfs(sg) isa Vector

show(isinteractive() ? stdout : devnull, sg)
show(isinteractive() ? stdout : devnull, MIME("text/plain"), sg)

@inferred sino_geom_help()

pg = sino_geom_plot_grids(plot)
ps = sino_geom_show()

plot(pg..., ps..., layout=(2,4))
prompt()
plot(pl[1:4]...)
