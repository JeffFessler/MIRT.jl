# kspace.jl

using MIRT: ir_mri_kspace_ga_radial

using Plots: plot, plot!, scatter
using Test: @test, @inferred


Nro = 30
delta_ro = 1
k = @inferred ir_mri_kspace_ga_radial(Nro=Nro, Nspoke=13, start=pi,
	shift=-0.5, delta_ro=delta_ro)
@test size(k) == (30, 13, 2)

kx = k[:,:,1]
ky = k[:,:,2]
scatter(kx, ky, label="", aspect_ratio=1,
	markerstrokecolor=:auto, markersize=3, markershape=:circle)
kmax = delta_ro * Nro / 2
plot!(xlim = [-1,1]*kmax)
plot!(ylim = [-1,1]*kmax)
plot!(xtick = (-1:1) * kmax, xlabel = "k_x")
plot!(ytick = (-1:1) * kmax, ylabel = "k_y")
