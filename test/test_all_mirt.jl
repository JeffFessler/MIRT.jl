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
)


"""
test_all_mirt()

run all MIRT tests
"""
function test_all_mirt()

	@testset "all MIRT" begin
		for root in srclist
		#	@show root
			@testset "$root" begin
				include("../src/" * root * "/z-test.jl")
			end
		end
	end

	true
end
