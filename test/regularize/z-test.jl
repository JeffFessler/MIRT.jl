# regularize/z-test.jl

using Test: @testset

list = [
"Aodwt.jl"
"diffl.jl"
"diffs.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
