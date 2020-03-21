#=
mri_trajectory.jl
based on mri_trajectory.m
Jing Dong
2020-03

using MIRT: image_geom_mri, ndgrid, mri_kspace_spiral, ir_mri_kspace_ga_radial, prompt
=#

export mri_trajectory

using Random: seed!
using Plots; default(markerstrokecolor=:auto, markersize=2)
using Test: @test


"""
    kspace, omega, wi = mri_trajectory( ; ktype, N, fov, arg_wi, kwargs...)

Generate kspace trajectory samples and density compensation functions.

option
* `ktype::Symbol` k-space trajectory type; default `:radial`
* `N::Dims` target image size; default (32,30)
* `fov` field of view in x and y (and z); default [250,250] mm
* `arg_wi` options to pass to `ir_mri_density_comp` - not yet done
* `kwargs` options for the specific trajectory

out
* `kspace [Nk 2|3]` kspace samples in units 1/fov
* `omega [Nk 2|3]` trajectory samples over [-π,π)
* `wi [Nk 1]` (optional) density compensation factors

trajectory types:
`:cartesian`
`:radial`
`:cart_y_2`
`:random`
`:half8`
`:epi_sin`
`:spiral0 :spiral1 :spiral3`
`:rosette3`
`:epi_under`
`:gads` (emulate golden-angle data sharing per winkelmann:07:aor)
"""
function mri_trajectory( ;
		ktype::Symbol = :radial,
		N::Dims = (32,30),
		fov::AbstractVector{<:Number} = 250 * ones(length(N)),
		wi::Union{AbstractArray{<:Real},Nothing} = nothing,
		arg_wi::Any = nothing, # not done
		kwargs...) # optional arguments for specific trajectories

	length(N) ∉ (2,3) && throw("only 2d and 3d done")

	# default density compensation, which works for ordinary Cartesian (only)
	wi = nothing

	# trajectory choices

	if (ktype == :cartesian)
		o1 = ((0:(N[1]-1))/N[1] .- 0.5)*2π
		o2 = ((0:(N[2]-1))/N[2] .- 0.5)*2π
		if (length(N) == 2)
			o1, o2 = ndgrid(o1, o2)
			omega = [o1[:] o2[:]]
		end

		if (length(N) == 3)
			o3 = ((0:(N[3]-1))/N[3] .- 0.5)*2π
			o1, o2, o3 = ndgrid(o1, o2, o3)
			omega = [o1[:] o2[:] o3[:]]
		end

		wi = 1 / prod(fov)
	end

	if (ktype == :epi_under)
		omega, wi = mri_trajectory_epi_under(N, fov ; kwargs...)
	end

	if (ktype == :gads)
		omega, wi = mri_trajectory_gads(N, fov ; kwargs...)
	end

	if (ktype == :radial)
@show N, fov
@show N[1:2], fov[1:2]
		omega, wi = mri_trajectory_radial(N[1:2], fov[1:2] ; kwargs...)

		if (length(N) == 3) # 'stack of radials'
			omega = mri_trajectory_stack(omega, N[3])
			wi = repeat(wi, N[3], 1)
			wi = wi[:] / fov[3]
		end
	end

	if (ktype == :rosette3)
		omega, wi = mri_trajectory_rosette3(N, fov ; kwargs...)
	end

	# echo-planar with sinusoid:
	if (ktype == :epi_sin)
		tmp = ( ; oversample=1) -> oversample
		oversample = tmp( ; kwargs...)
		Npt = round(Int, oversample * prod(N))
		t = (0:(Npt-1))/Npt
		omega = [π*sin.(2π*t*N[2]/2) t*2π .- π]
	end

	# bad spiral:
	if (ktype == :spiral0)
		Nspiral = 400
		omega = LinRange(0, 10*2π, Nspiral)
		omega = π * [cos.(omega) sin.(omega)] .* [omega omega]/maximum(omega)
		wi = abs.(omega[:,1] + 1im * omega[:,2]) # simple |r| weighting
	end

	# crude spiral:
	if (ktype == :spiral1)
		Nspiral = round(Int, prod(N) * π/4)
		omega = LinRange(0, N[1]*2π, Nspiral)
		max_omega = maximum(omega)
		omega = π * [cos.(omega) sin.(omega)].* [omega omega]/max_omega
		wi = abs.(omega[:,1] + 1im * omega[:,2]) # simple |r| weighting
	end

	# 3T spiral:
	if (ktype == :spiral3)
		kspace, omega = mri_kspace_spiral(N = maximum(N[1:2]), fov = maximum(fov[1:2]))

		if length(N) == 3 #stack of these spirals
			omega = mri_trajectory_stack(omega, N[3])
		end

		wi = abs.(omega[:,1] + 1im * omega[:,2]) # simple |r| weighting
		wi = wi / (fov[1]*fov[2]) # approximate scaling
		if length(N) == 3
			length(fov) != 3 && throw("fov not 3D")
			wi = wi / fov[3] # cartesian-style weighting in z
		end
	end

	# random
	if (ktype == :random)
		seed!(0)
		omega = (rand(prod(N), 2) .- 0.5) * 2π
	end

	# 2D half Cartesian + 8 rows
	if (ktype == :half8)
		extra = 8
		o1 = ((0:(N[1]-1))/N[1] .- 0.5)*2π
		o2 = (-N[2]/2:extra)/N[2] * 2π
		o1,o2 = ndgrid(o1, o2)
		omega = [o1[:] o2[:]]
	end

	# 2D FT, undersampled in "y" (phase encode) direction
	if (ktype == :cart_y_2)
		under = 2
		o1 = ((0:(N[1]/1-1))/(N[1]/1) .- 0.5) * 2π
		o2 = ((0:(N[2]/2-1))/(N[2]/under) .- 0.5) * 2π
		o1,o2 = ndgrid(o1, o2)
		omega = [o1[:] o2[:]]
	end

	# convert to physical units
	kspace = zeros(Float32, size(omega))
	for id = 1:length(N)
		dx = fov[id] / N[id]
		kspace[:,id] = omega[:,id] / (2π) / dx
	end

	# wi = ir_mri_density_comp(kspace, arg_wi{:}) # todo

	return kspace, omega, wi
end


# mri_trajectory_stack()
# make 3D "stack of ..." trajectory from a 2D trajectory
function mri_trajectory_stack(omega2, N3::Int)
	o3 = ((0:(N3-1))/N3 .- 0.5)*2π
	o3 = repeat(o3, 1, size(omega2,1)) # [N3,N12]
	o3 = o3'[:] # [N12*N3]
	omega = repeat(omega2, N3) # [N12*N3,2]
	omega = [omega o3] # [N12*N3,3]
	return omega
end


# mri_trajectory_epi_under()
# EPI, with optional under-sampling
function mri_trajectory_epi_under(N, fov ;
		under::Int = 2, # default is every other phase encode
		samp::AbstractVector{<:Bool} = rem.(0:(N[2]-1), under) .== 0,
	)
	nx = N[1]
	ny = N[2]
	omx = (-nx/2:nx/2-1) / nx * 2π # [-π,π) in x
	omy = (-ny/2:ny/2-1) / ny * 2π # [-π,π) in y
	omega = zeros(0,2)
	for iy in findall(samp)
		omega = [omega; [omx omy[iy]*ones(nx)]]
		omx = reverse(omx, dims = 1) # non-flyback
	end
	wi = ones(size(omega,1)) / prod(fov)
	return omega, wi
end


"""
    omega, wi = mri_trajectory_gads(N, fov ; ...)
emulate 2D golden angle radial sampling with data sharing

option:
`Nro` # of samples in each readout/spoke
`shift` shift along read-out due to gradient delays (stress)
`kmax_frac` fractions of maximum krad (0.5) for rings
`under` under-sampling factor for each annulus
"""
function mri_trajectory_gads(N, fov ;
		Nro::Int = maximum(N),
		delta_ro::Real = 1/Nro,
		shift::Real = -0.75,
		kmax_frac::Array = [0.20, 0.35, 0.501],
		nspoke::Array = floor.(Int, π * kmax_frac * Nro),
		under::Array = [1, 1, 0.6],
		start::AbstractArray = (0:2) * π/4
	)

	nspoke = floor.(Int, nspoke .* under)
	if true # make fibonacci for more uniform coverage
		phi = (1 + sqrt(5)) / 2
		n = round.(log.(nspoke * sqrt(5) .+ 0.5) / log(phi))
		nspoke = phi.^n / sqrt(5) .+ 0.5
		nspoke = floor.(Int, nspoke)
	end
	nring = length(nspoke)
	kmax_frac = [0; kmax_frac]
	omega = zeros(0,2)
	for ir=1:nring
		#todo: not sure here
		kspace = ir_mri_kspace_ga_radial(Nspoke = nspoke[ir],
			Nro = Nro, delta_ro = delta_ro, shift = shift, start = start[ir])

		kspace = reshape(kspace, :, 2)
		krad = sqrt.(sum(kspace.^2, dims = 2)[:])
		good = (kmax_frac[ir] .<= krad) .& (krad .< kmax_frac[ir+1])
		kspace = kspace[good, :]
		omega = [omega; 2π*kspace]
	end
	wi = [] # no default DCF
	return omega, wi
end


"""
    mri_trajectory_radial()

option:
* `na_nr` default ensures proper sampling at edge of k-space
* `na` angular spokes; default: na_nr * nr
* `nr` radial samples per spoke
* `ir` default: `0:nr`

todo: generalize to 3D using barger:02:trc
"""
function mri_trajectory_radial(N::Dims, fov ;
		nr::Int = round(Int, maximum(N)/2),
		ir::AbstractArray{Int} = 0:nr,
		na_nr::Real = 2π,
		na::Int = 4*ceil(Int, na_nr * nr/4), # mult of 4
	)

	om = ir/nr * π
	ang = (0:na-1)/na * 2π
	om, ang = ndgrid(om, ang) # [nr+1, na]
	omega = [(om .* cos.(ang))[:] (om .* sin.(ang))[:]]

	# density compensation factors based on "analytical" voronoi
	fov[2] != fov[1] && @warn("fov not square")
	du = 1/fov[1] # assume this radial sample spacing
	wi = π * du^2 / na * 2 * ir # see lauzon:96:eop, joseph:98:sei
	wi = collect(wi)
	wi[ir .== 0] .= π * (du/2)^2 / na

	wi = repeat(wi, na)[:]
	return omega, wi
end


"""
    mri_trajectory_rosette3(N, fov ; ...)
3d rosette, with default parameters from bucholz:08:miw

options:
omax: maximum omega
nt : time samples (65.536 ms for 4 usec dt)
dt : time sample spacing (4 usec)
ti : time samples
"""
function mri_trajectory_rosette3(N, fov ;
		f1::Real = 211,
		f2::Real = 117.13,
		f3::Real = length(N) == 3 ? 73.65 : 0,
		omax::Real = π,
		nshot::Int = 32, # todo: shots vs arms
		nt::Int = length(N) == 3 ? 16385 : 2^9, # time samples (65.536 ms for 4 usec dt)
		dt::Real = 4e-6, # time sample spacing (4 usec)
		ti::AbstractVector = (0:nt-1) * dt, # time samples
	)

	tmp = 2π * ti
	p1 = f1 * tmp
	p2 = f2 * tmp
	p3 = f3 * tmp
	kx = omax * sin.(p1) .* cos.(p2) .* cos.(p3)
	ky = omax * sin.(p1) .* sin.(p2) .* cos.(p3)
	kz = omax * sin.(p1) .* sin.(p3)
	omega = [kx ky kz]
	for is = 1:nshot-1 # n-shot, rotate kx,ky by 2π / N
		ang = is * 2π / nshot
		c = cos(ang)
		s = sin(ang)
		ox = c * kx + s * ky
		oy = -s * kx + c * ky
		omega = [omega; [ox oy kz]]
	end
	wi = omax^3 * abs.(sin.(p1).^2 .* cos.(p1) .* cos.(p3)) # bucholz:08:miw
	if length(N) == 2
		omega = omega[:,1:2]
	end
	return omega, wi
end


"""
    mri_trajectory(:test)
self test
"""
function mri_trajectory(test::Symbol)
	test != :test && throw(DomainError(test, "Not valid"))
	ig = image_geom_mri(nx = 2^5, ny = 2^5-4, fov = 250) # 250 mm FOV
	N = ig.dim

	ktype3 = [:cartesian, :spiral3, :radial] # 3D
	for ktype in ktype3 # 3D tests
		@show ktype
		kspace, omega, wi = mri_trajectory( ; ktype=ktype, N=(N...,4))
		@test size(kspace,2) == size(omega,2) == 3
		@test size(kspace,1) == size(omega,1)
		@test length(wi) ∈ (1, size(kspace,1))
		scatter(omega[:,1]/π, omega[:,2]/π, omega[:,3]/π, label="$ktype")
	end
	prompt()

	ktypes = [:spiral3, :spiral1, :spiral0, :radial, :gads, :rosette3,
		:epi_sin, :epi_under, :cartesian, :half8, :cart_y_2, :random,]
	# arg_wi = [:voronoi]
	plots = Any[]
	for ktype in ktypes
		@show ktype
		args =
			(ktype == :epi_sin) ? Dict(:oversample => 1/1,) :
			(ktype == :radial) ? Dict(:nr => 10, :na_nr => π/2) :
				Dict{Symbol,Nothing}()

		kspace, omega, wi = mri_trajectory( ;
			ktype=ktype, N=N, fov = ig.fovs, args...)
		p = scatter(omega[:,1]/π, omega[:,2]/π, aspect_ratio=1,
			xlim = [-1,1]*1.1, xtick = -1:1, xlabel = "omega1/π",
			ylim = [-1,1]*1.1, ytick = -1:1, ylabel = "omega2/π",
			title = string(ktype), label="")
		push!(plots, p)
	end
	plot(plots...); gui()

	return true
end
