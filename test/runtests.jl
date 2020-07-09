# test/runtests.jl

using MIRT
using Plots: default, plot
using Test: @testset, detect_ambiguities

if !isinteractive() # for usual noninteractive test, do not show plots
	default(show=false, reuse=true) # trick from Plots/runtests.jl
	ENV["GKSwstype"] = "100" # also from Plots/runtests.jl
end


"""
`function mirt_setup_plot()`

MIRT tests use Plots (with default GR backend).
When running on a remote linux server over X11 I get error messages from GKS.
The following kludge (?) tries to overcome this problem.
"""
function mirt_setup_plot()

	try
		plot(ones(3))
	catch
		if Sys.islinux
			ENV["GKS_WSTYPE"] = "x11"
		else
			@warn "Plot not working for Sys.KERNEL=$(Sys.KERNEL); test may fail"
		end
	end

	try
		plot(ones(3))
	catch
		@warn "Plot still not working; test may fail"
	end

	nothing
end

mirt_setup_plot()


# test_all_mirt
# Run all MIRT tests.
#To avoid plotting during tests, first do `default(show=false)`

# todo: old way (to cut):

srclist = (
"../data",
#"algorithm",
#"fbp",
#"io",
#"mri",
#"nufft",
"plot",
"regularize",
"system",
"utility",
)


basedir = dirname(pathof(MIRT))
#	@testset "all MIRT" begin
for root in srclist
	@testset "$root" begin
		tmp = joinpath(basedir, root, "z-test.jl")
	#	tmp = joinpath(root, "z-test.jl")
	#	printstyled(tmp * "\n", color=:blue, bold=true)
		include(tmp)
	end
end


# new way:
srclist = (
#"../data",
"algorithm",
"fbp",
"io",
"mri",
"nufft",
#"plot",
#"regularize",
#"system",
#"utility",
)

for root in srclist
	tmp = joinpath(root, "z-test.jl")
#	printstyled(tmp * "\n", color=:blue, bold=true)
	include(tmp)
end

@test length(detect_ambiguities(MIRT)) == 0
