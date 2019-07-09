# test/runtests.jl

using MIRT

using Plots: default

if !isinteractive()
	default(show=false, reuse=true) # trick from Plots/runtests.jl
end

# the following line is a trick that avoid plots during tests
# that take time and make CI output messages about "GKS"
# actually this line is no longer needed because of default() above!
# jim(:show, false)

mirt_setup_plot() # see src/plot/setup.jl

test_all_mirt()

# jim(:show, true) # return to default state where images are shown
