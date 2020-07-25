#=
sino_plot.jl
sinogram geometry plots
2020-07-24, Jeff Fessler, University of Michigan
=#

export sino_geom_plot!, sino_geom_plot_grids

# using MIRT: jim, image_geom, ImageGeom
# using MIRT: SinoGeom, sino_geom,
# using MIRT: SinoPar, SinoFan, SinoMoj, SinoFanArc, SinoFanFlat


"""
    sino_geom_plot_grid(sg, Plots.scatter) <- trick to avoid Plots.dependency
Scatter plot of `(r,phi)` sampling locations from `sg.grid`
"""
function sino_geom_plot_grid(sg::SinoGeom, scatter::Function)
	(r, phi) = sg.grid
	ylim = [min(0, rad2deg(minimum(phi))), max(360, rad2deg(maximum(phi)))]
	rmax = ceil(maximum(abs.(r))/10, digits=0)*10
	scatter(r, rad2deg.(phi), label="", markersize=1, markerstrokecolor=:auto,
		markershape=:circle, linewidth=0, # emulate scatter
		ylim = ylim, xlim = [-1,1]*rmax, # ylabel="ϕ",
		title="$(typeof(sg))", xtick=(-1:1)*rmax, ytick=[0,360])
end


"""
    sino_geom_plot_grids()
scatter plot of (r,phi) sampling locations for all geometries
"""
function sino_geom_plot_grids(
	plot::Function ;
	orbit::Real = 360, down::Int = 30,
)
	geoms = (
		sino_geom(:par ; nb = 888, na = 984, down=down, d = 0.5, orbit=orbit,
			offset = 0.25),
		sino_geom(:fan ; nb = 888, na = 984, d = 1.0, orbit = orbit,
			offset = 0.75, dsd = 949, dod = 408, down=down),
		sino_geom(:fan ; nb = 888, na = 984, d = 1.0, orbit = orbit,
			offset = 0.75, dsd = 949, dod = 408, down=down,
			dfs = Inf, source_offset = 0.7), # flat fan
		sino_geom(:moj ; nb = 888, na = 984, down=down, d = 1.0, orbit=orbit,
			offset = 0.25),
	)

	ngeom = length(geoms)
	pl = Array{Any}(undef, ngeom)

	for ii=1:ngeom
		sg = geoms[ii]
		pl[ii] = sg.plot_grid(plot)
	end
	return pl
end


"""
    sino_geom_plot!(sg, plot! ; ig)
Picture of the source position / detector geometry
"""
function sino_geom_plot!(
	sg::SinoGeom,
	plot!::Function ; # trick to avoid Plot dependency in MIRT
	ig::Union{Nothing,ImageGeom} = nothing,
)
	plot!(aspect_ratio=1)

	scat!(args... ; kwargs...) =
		plot!(args... ; kwargs..., markerstrokecolor=:auto, linewidth=0)

	xmax = sg.rfov; xmin = -xmax; (ymin,ymax) = (xmin,xmax)
	if !isnothing(ig)
#		plot!(jim(ig.x, ig.y, ig.mask[:,:,1], clim=(0,1))) # todo: jim!
		xmin = minimum(ig.x); xmax = maximum(ig.x)
		ymin = minimum(ig.y); ymax = maximum(ig.y)
	end
	plot!([xmax, xmin, xmin, xmax, xmax],
		[ymax, ymax, ymin, ymin, ymax], color=:green, label="")
	plot!(xtick=round.([xmin, 0, xmax], digits=0))
	plot!(ytick=round.([ymin, 0, ymax], digits=2))

	θ = LinRange(0, 2*pi, 1001)
	rfov = sg.rfov
	scat!([0], [0], marker=:circle, label="")
	plot!(rfov * cos.(θ), rfov * sin.(θ), color=:magenta, label="") # rfov circle
	rfov = round(sg.rfov, digits=1)
	plot!(xlabel="x", ylabel="y", title = "$(typeof(sg))")

#=
	if sg isa SinoPar
	end
=#

	if sg isa SinoFan
		x0 = 0
		y0 = sg.dso
		t = LinRange(0, 2π, 100)
		rot = sg.ar[1]
		rot = [cos(rot) -sin(rot); sin(rot) cos(rot)]
		p0 = rot * [x0; y0]
		pd = rot * [sg.xds'; sg.yds'] # detector points

		tmp = sg.ar .+ π/2 # trick: angle beta defined ccw from y axis
		scat!([p0[1]], [p0[2]], color=:blue, marker=:square, label="") # source
		plot!(sg.dso * cos.(t), sg.dso * sin.(t), color=:cyan, label="") # source circle
		scat!(sg.dso * cos.(tmp), sg.dso * sin.(tmp),
			color=:blue, marker=:circle, markersize=2, label="") # source points
		scat!(vec(pd[1,:]), vec(pd[2,:]), marker=:circle,
			color=:orange, markersize=1, label="") # detectors

		plot!([pd[1,1], p0[1], pd[1,end]], [pd[2,1], p0[2], pd[2,end]],
			color=:red, label="")
		plot!(title="$(typeof(sg))")
	end

	if sg isa SinoMoj
		θ = LinRange(0, 2π, 100)
		rphi = sg.nb/2 * sg.d_moj.(θ)
		plot!(rphi .* cos.(θ), rphi .* sin.(θ), color=:blue, label="")
	#	rmax = maximum(sg.s)
	#	axis([-1 1 -1 1] * max([rmax ig.fov/2]) * 1.1)
	end

	plot!()
end
