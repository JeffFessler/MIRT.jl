# test/test_all_mirt.jl

using Test

srclist = (
"../data",
"algorithm",
"fbp",
"io",
"plot",
"regularize",
"system",
"utility",
)


"""
test_all_mirt()

run all MIRT tests
"""
function test_all_mirt()

	basedir = dirname(pathof(MIRT))
	@testset "all MIRT" begin
		for root in srclist
		#	@show root
			@testset "$root" begin
				tmp = joinpath(basedir, root, "z-test.jl")
			#	printstyled(tmp * "\n", color=:blue, bold=true)
				include(tmp)
			end
		end
	end

	true
end
