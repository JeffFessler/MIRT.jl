#=
sino_geom.jl
sinogram geometry for 2D tomographic image reconstruction
2019-07-01, Jeff Fessler, University of Michigan
=#

export MIRT_sino_geom, sino_geom
export sino_geom_help, sino_geom_plot_grids, sino_geom_show, sino_geom_test

# using MIRT: jim, image_geom, MIRT_image_geom, prompt
using Plots: Plot, plot!, plot, scatter!, gui
using Test: @test, @test_throws, @inferred


"""
    MIRT_sino_geom
struct to describe a 2D sinogram geometry
"""
struct MIRT_sino_geom
	how::Symbol				# :par | :moj | :fan 
	units::Symbol			# :nothing | :mm | :cm etc.
	nb::Int					# # of "radial" samples, aka ns
	na::Int					# # of angular samples
	d::Float32				# aka dr or ds, "radial" sample spacing
							# (is dx for mojette, pixels must be square)
	orbit::Float32			# [degrees]
	orbit_start::Float32	# [degrees]
	offset::Float32			# sample offset, cf offset_r or offset_s [unitless]
	strip_width::Float32	# 

	# for fan:
	source_offset::Float32	# same units as d, etc., e.g., [mm]
							# use with caution!
	dsd::Float32			# dis_src_det, Inf for parallel beam
	dod::Float32			# dis_iso_det
#	dso::Float32			# dis_src_iso = dsd-dod, Inf for parallel beam
	dfs::Float32			# distance from focal spot to source
end


"""
    sino_geom_help()
"""
function sino_geom_help( ; io::IO = isinteractive() ? stdout : devnull )
	print(io, "$(basename(@__FILE__)) propertynames:\n\t")
	print(io, Tuple(sort([propertynames(sino_geom(:par))...])))

	print(io,
	"\n
	Derived values

	sg.dim			dimensions: (nb,na)
	sg.ds|dr		radial sample spacing (NaN for :moj)
	sg.s			[nb] s sample locations
	sg.w			(nb-1)/2 + offset ('middle' sample position)
	sg.ad			source angles [degrees]
	sg.ar			source angles [radians]
	sg.ones			ones(Float32, nb,na)
	sg.zeros		zeros(Float32, nb,na)
	sg.rfov			radial fov
	sg.xds			[nb] center of detector elements (beta=0)
	sg.yds			[nb] ''
	sg.grid			(rg, phigrid) [nb na] parallel-beam coordinates
	sg.plot_grid	plot sg.grid

	For mojette:

	sg.d_ang		[na]

	For fan beam:

	sg.gamma		[nb] gamma sample values [radians]
	sg.gamma_max	half of fan angle [radians]
	sg.dso			# dsd - dod, Inf for parallel beam

	Methods

	sg.down(down)		reduce sampling by integer factor
	sg.shape(sino)		reshape sinograms into array [nb na :]
	sg.unitv(;ib,ia)	unit 'vector' with single nonzero element
	sg.taufun(x,y)		projected s/ds for each (x,y) pair [numel(x) na]
	sg.plot(;ig)		plot system geometry (most useful for fan)
	\n")
end


"""
    function sg = sino_geom(...)

Constructor for `MIRT_sino_geom`

Create the "sinogram geometry" structure that describes the sampling
characteristics of a given sinogram for a 2D parallel or fan-beam system.
Using this structure facilitates "object oriented" code.
(Use `ct_geom()` instead for 3D axial or helical cone-beam CT.)

in
- `how::Symbol`	`:fan` (fan-beam) | `:par` (parallel-beam) | `:moj` (mojette)

options for all geometries (including parallel-beam):
- `units::Symbol`	e.g. `:cm` or `:mm`; default: :none
- `orbit_start`		default: 0
- `orbit`			[degrees] default: `180` for parallel / mojette
					and `360` for fan
					can be `:short` for fan-beam short scan
- `down::Int`		down-sampling factor, for testing

- `nb`				# radial samples cf `nr` (i.e., `ns` for `:fan`)
- `na`				# angular samples (cf `nbeta` for `:fan`)
- `d`				radial sample spacing; cf `dr` or `ds`; default 1
					for mojette this is actually `dx`
- `offset`			cf `offset_r` `channel_offset` unitless; default 0
			(relative to centerline between two central channels).
			Use 0.25 or 1.25 for "quarter-detector offset"
- `strip_width`		detector width; default: `d`

options for fan-beam
- `source_offset`		same units as d; use with caution! default 0
fan beam distances:
- `dsd`		cf `dis_src_det`	default: `Inf` (parallel beam)
- `dod`		cf `dis_iso_det`	default: `0`
- `dfs`		cf `dis_foc_src`	default: `0` (3rd generation CT arc),
				use `Inf` for flat detector

out
- `sg::MIRT_sino_geom`	initialized structure

See also
- `sino_geom_help()` help on methods
- `sino_geom_plot_grids()` show sampling
- `sino_geom_show()` show the geometries
- `sino_geom_test()` self test

Jeff Fessler, University of Michigan
"""
function sino_geom(how::Symbol; kwarg...)
	if how === :par
		sg = sino_geom_par( ; kwarg...)
	elseif how === :fan
		sg = sino_geom_fan( ; kwarg...)
	elseif how === :moj
		sg = sino_geom_moj( ; kwarg...)
	elseif how === :ge1
		sg = sino_geom_ge1( ; kwarg...)
#=
	elseif how === :hd1
		sg = sino_geom_hd1( ; kwarg...)
	elseif how === :revo1fan
		tmp = ir_fan_geom_revo1(type)
		sg = sino_geom(:fan, tmp{:}, varargin{:})
=#
	else
		throw("unknown sino type $how")
	end

	return sg
end


"""
    sg = downsample(sg, down)
down-sample (for testing with small problems)
"""
function downsample(sg::MIRT_sino_geom, down::Int)
	down == 1 && return sg

	nb = 2 * round(Int, sg.nb / down / 2) # keep it even
	na = round(Int, sg.na / down)

	return MIRT_sino_geom(sg.how, sg.units,
		nb, na, sg.d * down, sg.orbit, sg.orbit_start, sg.offset,
		sg.strip_width * down,
		sg.source_offset, sg.dsd, sg.dod, sg.dfs)
end


"""
    sg = sino_geom_over(sg, over::Int)
over-sample in "radial" dimension
Probably not meaningful for mojette sampling because d=dx.
"""
function sino_geom_over(sg::MIRT_sino_geom, over::Int)
	if over == 1
		return sg
	end

	return MIRT_sino_geom(sg.how, sg.units,
		sg.nb * over, sg.na, sg.d / over,
		sg.orbit, sg.orbit_start, sg.offset * over,
		sg.strip_width / over,
		sg.source_offset, sg.dsd, sg.dod, sg.dfs)
end


"""
    sg = sino_geom_fan()
"""
function sino_geom_fan( ;
		units::Symbol = :none,
		nb::Int = 128,
		na::Int = 2 * floor(Int, nb * pi/2 / 2),
		d::Real = 1,
		orbit::Union{Symbol,Real} = 360, # [degrees]
		orbit_start::Real = 0,
		strip_width::Real = d,
		offset::Real = 0,
		source_offset::Real = 0,
		dsd::Real = 4*nb*d,	# dis_src_det
	#	dso::Real = [],		# dis_src_iso
		dod::Real = nb*d,	# dis_iso_det
		dfs::Real = 0,		# dis_foc_src (3rd gen CT)
		down::Int = 1,
	)

	dfs != 0 && !isinf(dfs) && throw("dfs $dfs") # must be 0 or Inf

	if orbit === :short # trick
		sg_tmp = MIRT_sino_geom(:fan, units,
			nb, na, d, 0, orbit_start, offset, strip_width,
			source_offset, dsd, dod, dfs)
		orbit = sg_tmp.orbit_short
	end
	isa(orbit, Symbol) && throw("orbit :orbit")

	sg = MIRT_sino_geom(:fan, units,
		nb, na, d, orbit, orbit_start, offset, strip_width,
		source_offset, dsd, dod, dfs)

	return downsample(sg, down)
end


"""
    sg = sino_geom_par( ... )
"""
function sino_geom_par( ;
		units::Symbol = :none,
		nb::Int = 128,
		na::Int = 2 * floor(Int, nb * pi/2 / 2),
		down::Int = 1,
		d::Real = 1,
		orbit::Real = 180, # [degrees]
		orbit_start::Real = 0,
		strip_width::Real = d,
		offset::Real = 0,
	)

	sg = MIRT_sino_geom(:par, units,
		nb, na, d, orbit, orbit_start, offset, strip_width,
		0, 0, 0, 0)

	return downsample(sg, down)
end


"""
    sg = sino_geom_moj( ... )
"""
function sino_geom_moj( ;
		units::Symbol = :none,
		nb::Int = 128,
		na::Int = 2 * floor(Int, nb * pi/2 / 2),
		down::Int = 1,
		d::Real = 1, # means dx for :moj
		orbit::Real = 180, # [degrees]
		orbit_start::Real = 0,
		strip_width::Real = d, # ignored ?
		offset::Real = 0,
	)

	sg = MIRT_sino_geom(:moj, units,
		nb, na, d, orbit, orbit_start, offset, strip_width,
		0, 0, 0, 0)

	return downsample(sg, down)
end


"gamma for general finite dfs (rarely used)"
function sino_geom_gamma_dfs(sg)
	dis_foc_det = sg.dfs + sg.dsd
	alf = sg.s / dis_foc_det
	atan.(dis_foc_det * sin.(alf), dis_foc_det * cos.(alf) .- sg.dfs)
	# equivalent to s/dsd when dfs=0
end


"""
    sino_geom_gamma()
gamma sample values for :fan
"""
function sino_geom_gamma(sg)
	return	sg.dfs == 0 ? sg.s / sg.dsd : # 3rd gen: equiangular
			isinf(sg.dfs) ? atan.(sg.s / sg.dsd) : # flat
			sino_geom_gamma_dfs(sg) # general
end


"""
    sino_geom_rfov()
radial FOV
"""
function sino_geom_rfov(sg)
	return	sg.how === :par ? maximum(abs.(sg.r)) :
			sg.how === :fan ? sg.dso * sin(sg.gamma_max) :
			sg.how === :moj ? sg.nb/2 * minimum(sg.d_ang) : # todo: check
				throw("bad how $(sg.how)")
end


"""
    sino_geom_taufun()
projected `s/ds`, useful for footprint center and support
"""
function sino_geom_taufun(sg, x, y)
	size(x) != size(y) && throw("bad x,y size")
	x = vec(x)
	y = vec(y)
	if sg.how === :par || sg.how === :moj # todo: check
		ar = sg.ar' # row vector, for outer-product
		tau = (x * cos.(ar) + y * sin.(ar)) / sg.dr
	elseif sg.how === :fan
		b = sg.ar' # row vector, for outer-product
		xb = x * cos.(b) + y * sin.(b)
		yb = -x * sin.(b) + y * cos.(b)
		tangam = (xb .- sg.source_offset) ./ (sg.dso .- yb) # e,tomo,fan,L,gam
		if sg.dfs == 0 # arc
			tau = sg.dsd / sg.ds * atan.(tangam)
		elseif isinf(sg.dfs) # flat
			tau = sg.dsd / sg.ds * tangam
#		else
#			throw("bad dfs $(sg.dfs)")
		end
#	else
#		throw("bad how $(sg.how)")
	end
	return tau
end


"""
    sino_geom_xds()
center positions of detectors (for beta = 0)
"""
function sino_geom_xds(sg)
	if sg.how === :par
		xds = sg.s
	elseif sg.how === :moj
		xds = sg.s # todo: really should be angle dependent
	elseif sg.how === :fan
		if sg.dfs == 0 # arc
			gam = sg.gamma
			xds = sg.dsd * sin.(gam)
		elseif isinf(sg.dfs) # flat
			xds = sg.s
	#	else
	#		throw("bad dfs $(sg.dfs))")
		end
#	else
#		throw("bad how $(sg.how)")
	end
	return xds .+ sg.source_offset
end


"""
    sino_geom_yds()
center positions of detectors (for beta = 0)
"""
function sino_geom_yds(sg)

	if sg.how === :par
		yds = zeros(Float32, sg.nb)
	elseif sg.how === :moj
		yds = zeros(Float32, sg.nb)
	elseif sg.how === :fan
		if sg.dfs == 0 # arc
			gam = sg.gamma
			yds = sg.dso .- sg.dsd * cos.(gam)
		elseif isinf(sg.dfs) # flat
			yds = fill(-sg.dod, sg.nb)
	#	else
	#		throw("bad dfs $(sg.dfs))")
		end
#	else
#		throw("bad how $(sg.how)")
	end
	return yds
end


"""
    sino_geom_unitv()
sinogram with a single ray
"""
function sino_geom_unitv(sg::MIRT_sino_geom;
		ib=round(Int, sg.nb/2+1),
		ia=round(Int, sg.na/2+1))
	out = sg.zeros
	out[ib,ia] = 1
	return out
end


"""
    (rg, ϕg) = sino_geom_grid(sg::MIRT_sino_geom)

Return grids `rg` and `ϕg` (in radians) of size `[nb na]`
of equivalent *parallel-beam* `(r,ϕ)` (radial, angular) sampling positions,
for any sinogram geometry.
For parallel beam this is just `ndgrid(sg.r, sg.ar)`
but for fan beam and mojette this involves more complicated computations.
"""
function sino_geom_grid(sg::MIRT_sino_geom)

	if sg.how === :par
		return ndgrid(sg.r, sg.ar)

	elseif sg.how === :fan
		gamma = sg.gamma
		rad = sg.dso * sin.(gamma) + sg.source_offset * cos.(gamma)
		rg = repeat(rad, 1, sg.na) # [nb na]
		return (rg, gamma .+ sg.ar') # [nb na] phi = gamma + beta
	end

	# otherwise :moj (mojette)
	phi = sg.ar
	# trick: ray_spacing aka ds comes from dx which is sg.d for mojette
	wb = (sg.nb - 1)/2 + sg.offset
	dt = sg.d_ang # [na]
	pos = ((0:(sg.nb-1)) .- wb) * dt' # [nb na]
	return (pos, repeat(phi', sg.nb, 1))
end


"""
    show(io::IO, sg::MIRT_sino_geom)
    show(io::IO, ::MIME"text/plain", sg::MIRT_sino_geom)
"""
Base.show(io::IO, sg::MIRT_sino_geom) =
    print(io, "MIRT_sino_geom: $(sg.dim)")
function Base.show(io::IO, ::MIME"text/plain", sg::MIRT_sino_geom)
    ir_dump(io, sg)
end


# Extended properties

sino_geom_fun0 = Dict([
	(:help, sg -> sino_geom_help()),

	(:dim, sg -> (sg.nb, sg.na)),
	(:w, sg -> (sg.nb-1)/2 + sg.offset),
	(:ones, sg -> ones(Float32, sg.dim)),
	(:zeros, sg -> zeros(Float32, sg.dim)),

	(:dr, sg -> sg.how === :moj ? NaN : sg.d),
	(:ds, sg -> sg.how === :moj ? NaN : sg.d),
	(:r, sg -> sg.d * ((0:sg.nb-1) .- sg.w)),
	(:s, sg -> sg.r), # sample locations ('radial')

	(:gamma, sg -> sino_geom_gamma(sg)),
	(:gamma_max, sg -> maximum(abs.(sg.gamma))),
	(:orbit_short, sg -> 180 + 2 * rad2deg(sg.gamma_max)),
	(:ad, sg -> (0:sg.na-1)/sg.na * sg.orbit .+ sg.orbit_start),
	(:ar, sg -> deg2rad.(sg.ad)),

	(:rfov, sg -> sino_geom_rfov(sg)),
	(:xds, sg -> sino_geom_xds(sg)),
	(:yds, sg -> sino_geom_yds(sg)),
	(:dso, sg -> sg.dsd - sg.dod),
	(:grid, sg -> sino_geom_grid(sg)),
	(:plot_grid, sg -> sino_geom_plot_grid(sg)),

	# angular dependent d for :moj
	(:d_ang, sg -> sg.d * max.(abs.(cos.(sg.ar)), abs.(sin.(sg.ar)))),

	(:shape, sg -> ((x::AbstractArray{<:Number} -> reshape(x, sg.dim..., :)))),
	(:taufun, sg -> ((x,y) -> sino_geom_taufun(sg,x,y))),
	(:unitv, sg -> ((;kwarg...) -> sino_geom_unitv(sg; kwarg...))),
	(:plot, sg -> ((;ig=nothing) -> sino_geom_plot(sg, ig=ig))),

	# functions that return new geometry:

	(:down, sg -> (down::Int -> downsample(sg, down))),
	(:over, sg -> (over::Int -> sino_geom_over(sg, over))),

	])


# Tricky overloading here!

Base.getproperty(sg::MIRT_sino_geom, s::Symbol) =
		haskey(sino_geom_fun0, s) ? sino_geom_fun0[s](sg) :
		getfield(sg, s)

Base.propertynames(sg::MIRT_sino_geom) =
	(fieldnames(typeof(sg))..., keys(sino_geom_fun0)...)


"""
    sino_geom_plot_grid()
scatter plot of (r,phi) sampling locations from `sg.grid`
"""
function sino_geom_plot_grid(sg::MIRT_sino_geom)
	(r, phi) = sg.grid
	dfs = sg.how === :fan ? " dfs=$(sg.dfs)" : ""
	ylim = [min(0, rad2deg(minimum(phi))), max(360, rad2deg(maximum(phi)))]
	rmax = ceil(maximum(abs.(r))/10, digits=0)*10
	scatter(r, rad2deg.(phi), label="", markersize=1, markerstrokecolor=:auto,
		ylim = ylim, xlim = [-1,1]*rmax, # ylabel="ϕ",
		title="$(sg.how)$dfs", xtick=(-1:1)*rmax, ytick=[0,360])
end


"""
    sino_geom_plot_grids()
scatter plot of (r,phi) sampling locations for all geometries
"""
function sino_geom_plot_grids( ; orbit::Real = 360, down::Int = 30)
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
	pl = Array{Plot}(undef, ngeom)

	for ii=1:ngeom
		sg = geoms[ii]
		pl[ii] = sg.plot_grid
	end
	return pl
end


"""
    sino_geom_plot()
picture of the source position / detector geometry
"""
function sino_geom_plot(sg ; ig::Union{Nothing,MIRT_image_geom}=nothing)
	plot(aspect_ratio=1)

	xmax = sg.rfov; xmin = -xmax; (ymin,ymax) = (xmin,xmax)
	if !isnothing(ig)
		plot!(jim(ig.x, ig.y, ig.mask[:,:,1], clim=(0,1)))
		xmin = minimum(ig.x); xmax = maximum(ig.x)
		ymin = minimum(ig.y); ymax = maximum(ig.y)
	end
	plot!([xmax, xmin, xmin, xmax, xmax],
		[ymax, ymax, ymin, ymin, ymax], color=:green, label="")
	plot!(xtick=round.([xmin, 0, xmax], digits=0))
	plot!(ytick=round.([ymin, 0, ymax], digits=2))

	t = LinRange(0, 2*pi, 1001)
	rfov = sg.rfov
	scatter!([0], [0], marker=:circle, label="")
	plot!(rfov * cos.(t), rfov * sin.(t), color=:magenta, label="") # rfov circle
	rfov = round(sg.rfov, digits=1)
	plot!(xlabel="x", ylabel="y", title = "$(sg.how): rfov = $rfov")

#=
	if sg.how === :par
	end
=#

	if sg.how === :fan
		x0 = 0
		y0 = sg.dso
		t = LinRange(0, 2*pi, 100)
		rot = sg.ar[1]
		rot = [cos(rot) -sin(rot); sin(rot) cos(rot)]
		p0 = rot * [x0; y0]
		pd = rot * [sg.xds'; sg.yds'] # detector points

		tmp = sg.ar .+ pi/2 # trick: angle beta defined ccw from y axis
		scatter!([p0[1]], [p0[2]], color=:yellow, label="") # source
		plot!(sg.dso * cos.(t), sg.dso * sin.(t), color=:cyan, label="") # source circle
		scatter!(sg.dso * cos.(tmp), sg.dso * sin.(tmp), markerstrokecolor=:auto,
			color=:cyan, markersize=1, label="") # source points
		scatter!(vec(pd[1,:]), vec(pd[2,:]), markerstrokecolor=:auto,
			color=:yellow, markersize=1, label="") # detectors

		plot!([pd[1,1], p0[1], pd[1,end]], [pd[2,1], p0[2], pd[2,end]],
			color=:red, label="")
		plot!(title="$(sg.how): dfs = $(sg.dfs)")
	end

	if sg.how === :moj && false
		t = LinRange(0, 2*pi, 100)
		rmax = maximum(sg.s)
		rphi = sg.nb/2 * sg.d ./ (max(abs.(cos.(t)), abs.(sin.(t))))
		plot!(rphi .* cos.(t), rphi .* sin.(t), label="") # fov circle
	#	axis([-1 1 -1 1] * max([rmax ig.fov/2]) * 1.1)
	end

	plot!()
end


"""
    sino_geom_ge1()
sinogram geometry for GE lightspeed system
These numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl
"""
function sino_geom_ge1( ;
		na::Int = 984,
		nb::Int = 888,
		orbit::Union{Symbol,Real} = 360,
		units::Symbol = :mm, # default units is mm
		kwarg...)

	if orbit === :short
		na = 642 # trick: reduce na for short scans
		orbit = na/984*360
	end

	scale = units === :mm ? 1 :
			units === :cm ? 10 :
			throw("units $units")

	return sino_geom(:fan ; units=units,
			nb=nb, na=na,
			d = 1.0239/scale, offset = 1.25,
			dsd = 949.075/scale,
			dod = 408.075/scale,
			dfs = 0, kwarg...)
end
