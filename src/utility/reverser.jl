#=
reverser.jl
based on
https://stackoverflow.com/questions/27411401/julia-reverse-n-dimensional-arrays
=#

export reverser

using Test: @test

"""
    y = reverser(x, dims)
reverse array along specified dimensions (or all if unspecified)
"""
function reverser(x::AbstractArray, dims::AbstractVector{<:Int})
    y = copy(x)
    for d in dims
        y = reverse(y, dims=d)
    end
    return y
end

reverser(x::AbstractArray) = reverser(x, 1:ndims(x)) # all dimensions
reverser(x::AbstractArray, d::Int) = reverser(x, [d])


"""
    reverser(:test)
self test
"""
function reverser(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	@test reverser(1:3) == 3:-1:1
	@test reverser(1:3, 1) == 3:-1:1
	@test reverser((1:3)', 1) == (1:3)'
	@test reverser((1:3)', 2) == (3:-1:1)'
	true
end
