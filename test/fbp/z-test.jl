# fbp/z-test.jl

using Test: @testset

list = [
"rotate2d.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
