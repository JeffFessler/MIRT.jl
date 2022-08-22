# fbp/z-test.jl

using Test: @testset

list = [
"rotate2d.jl"
"sino_geom.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
