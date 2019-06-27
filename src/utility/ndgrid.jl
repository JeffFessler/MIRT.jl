#=
ndgrid.jl

I know that using ndgrid is not the "julian way" but when porting code
from matlab sometimes it is easier just to use it instead of refactoring
to do some other way.

2019-06-18, Jeff Fessler, University of Michigan
=#

using Test: @test


"""
`(xx,yy) = ndgrid(x::AbstractVector{<:Any}, y::AbstractVector{<:Any})`
todo - improve?
"""
function ndgrid(
		x::AbstractVector{<:Any},
		y::AbstractVector{<:Any})

	tmp = Iterators.product(x, y)
	return [p[1] for p in tmp], [p[2] for p in tmp]
end


"""
`(xx,yy,zz) = ndgrid(x::AbstractVector{<:Any}, y::..., z::...)`
todo - improve?
"""
function ndgrid(
		x::AbstractVector{<:Any},
		y::AbstractVector{<:Any},
		z::AbstractVector{<:Any})

	tmp = Iterators.product(x, y, z)
	return [p[1] for p in tmp], [p[2] for p in tmp], [p[3] for p in tmp]
end


"""
`ndgrid(:test)`
self test
"""
function ndgrid(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))

	@test ndgrid(1:3, 2:4) == ([1 1 1; 2 2 2; 3 3 3], [2 3 4; 2 3 4; 2 3 4])
	@test size(ndgrid(1:3, 4:5, 6:9)[3]) == (3,2,4)
	true
end
