# algorithm/z-test.jl

using Test: @testset

list = [
"general/pogm_restart.jl"
"general/ncg.jl"
"general/ogm_ls.jl"
"general/poweriter.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
