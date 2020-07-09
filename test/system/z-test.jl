# system/z-test.jl

using Test: @testset

list = [
"Afft.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
