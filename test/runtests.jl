# test/runtests.jl

using MIRT

# the following line is a trick that avoid plots during tests
# that take time and make CI output messages about "GKS"
jim(:show, false)

mirt_setup_plot() # see src/plot/setup.jl

test_all_mirt()

jim(:show, true) # return to default state where images are shown
