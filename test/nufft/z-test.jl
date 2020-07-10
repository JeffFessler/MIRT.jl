# nufft/z-test.jl

using Test: @testset

# nufft/z-list.jl

list = [
"dtft.jl"
"nufft.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
