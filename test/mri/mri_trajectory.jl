# mri_trajectory.jl


using MIRT: mri_trajectory, image_geom_mri
using MIRTjim: prompt

using Plots
using Test: @test


ig = image_geom_mri(nx = 2^5, ny = 2^5-4, fov = 250) # 250 mm FOV
N = ig.dim

ktype3 = [:cartesian, :spiral3, :radial] # 3D
for ktype in ktype3 # 3D tests
	isinteractive() && (@show ktype)
	kspace, omega, wi = mri_trajectory( ; ktype=ktype, N=(N...,4))
	@test size(kspace,2) == size(omega,2) == 3
	@test size(kspace,1) == size(omega,1)
	@test length(wi) ∈ (1, size(kspace,1))
	scatter(omega[:,1]/π, omega[:,2]/π, omega[:,3]/π, label="$ktype",
		markerstrokecolor=:auto, markersize=2,
	)
end
prompt()

ktypes = [:spiral3, :spiral1, :spiral0, :radial, :gads, :rosette3,
	:epi_sin, :epi_under, :cartesian, :half8, :cart_y_2, :random,]
# arg_wi = [:voronoi]
plots = Array{Any}(undef, length(ktypes))
for (i,ktype) in enumerate(ktypes)
	isinteractive() && (@show ktype)
	args =
		(ktype == :epi_sin) ? Dict(:oversample => 1/1,) :
		(ktype == :radial) ? Dict(:nr => 10, :na_nr => π/2) :
			Dict{Symbol,Nothing}()

	kspace, omega, wi = mri_trajectory( ;
		ktype=ktype, N=N, fov = ig.fovs, args...)
	plots[i] = scatter(omega[:,1]/π, omega[:,2]/π, aspect_ratio=1,
		xlim = [-1,1]*1.1, xtick = -1:1, xlabel = "omega1/π",
		ylim = [-1,1]*1.1, ytick = -1:1, ylabel = "omega2/π",
		title = string(ktype), label="",
		markerstrokecolor=:auto, markersize=1,
	)
end
plot(plots...); gui()
