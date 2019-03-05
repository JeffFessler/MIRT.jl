# test/runtests.jl

using Test
using MIRT

srclist = (
"../data",
"algorithm",
"fbp",
"io",
"plot",
"regularize",
"system",
)

@testset "all MIRT" begin
	for root in srclist
	@show root
		@testset "$root" begin
			include("../src/" * root * "/z-test.jl")
		end
	end
end
