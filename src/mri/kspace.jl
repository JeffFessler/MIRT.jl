# kspace.jl
# kspace sampling patterns
# 2018-06-13 Jeff Fessler, University of Michigan

using Plots: plot, plot!, scatter, gui
using Test


"""
`kspace = ir_mri_kspace_ga_radial(; Nro=?, Nspoke=?, ...)`

Generate k-space sampling pattern for "golden angle" radial sampling.

option
* `Nro:Int`		number of samples in each readout/spoke, default 256
* `Nspoke::Int`	number of spokes, default 1
* `start::Real`	first angle in series [radians], default pi/2
* `angle::Real`	angular spacing [radians], default GA
* `delta_ro::Real`	readout spacing, default 1/Nro
* `shift::Real`		shift due to gradient delays, default 0
		radial sample locations are `ir * delta_ro`
		where `ir = [-(Nro/2 - 1):1:Nro/2] + shift`
* `show::Bool`		plot k-space locations, default `false`

out
* `kspace`	`[Nro Nspoke 2]` (Float32)
			kx and ky k-space locations for Nspoke*Nro samples
			in interval (-0.5 0.5] for default shift, delta_ro
			so default units are "cycles / sample"

2015-07 Mai Le, original Matlab version
2015-07-04 Jeff Fessler, minor changes to Matlab version
"""
function ir_mri_kspace_ga_radial(;
		Nro::Int = 256,
		Nspoke::Int = 1,
		delta_ro::Real = 1/Nro,
		shift::Real = 0,
		angle::Real = pi*(sqrt(5)-1)/2, # golden angle [radians]
		start::Real = pi/2,
		show::Bool = false,
	)

	isodd(Nro) && throw(ArgumentError("Number along read out not even!"))

	rho = (-(Nro/2 - 1):1:Nro/2) * delta_ro # radial samples
	rho = rho .+ shift * delta_ro # apply gradient delay shift

	phi = (0:Nspoke-1) * angle .+ start # angles

	rho = Float32.(rho)
	phi = Float32.(phi)

	kx = rho * cos.(phi)' # [Nro Nspoke]
	ky = rho * sin.(-phi)'
	kspace = cat(dims=3, kx, ky) # [Nro Nspoke 2]

	if show
		scatter(kx, ky, label="", aspect_ratio=1)
		kmax = delta_ro * Nro / 2
		plot!(xlim = [-1,1]*kmax)
		plot!(ylim = [-1,1]*kmax)
		plot!(xtick = (-1:1) * kmax, xlabel = "k_x")
		plot!(ytick = (-1:1) * kmax, ylabel = "k_y")
		gui()
	end

	return kspace
end


"""
`ir_mri_kspace_ga_radial`
self test
"""
function ir_mri_kspace_ga_radial(test::Symbol)
	test != :test && throw(":test")
	k = ir_mri_kspace_ga_radial(Nro=30, Nspoke=13, start=pi,
		show=true, shift=-0.5, delta_ro=1)
	@test size(k) == (30, 13, 2)
	true
end


# ir_mri_kspace_ga_radial(:test)
