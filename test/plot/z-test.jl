# plot/z-test.jl

using Test: @testset

list = [
"jim.jl"
#"setup.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
