# test/runtests.jl

using MIRT

using Plots: default

if !isinteractive() # for usual noninteractive test, do not show plots
	default(show=false, reuse=true) # trick from Plots/runtests.jl
	ENV["GKSwstype"] = "100" # also from Plots/runtests.jl
end

mirt_setup_plot() # see src/plot/setup.jl

test_all_mirt()
