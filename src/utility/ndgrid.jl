#=
ndgrid.jl

I know that using ndgrid is not the "julian way" but when porting code
from matlab sometimes it is easier just to use it instead of refactoring
to do some other way.

The better way would be to use CartesianIndices
https://julialang.org/blog/2016/02/iteration

2019-06-18, Jeff Fessler, University of Michigan
=#

export ndgrid


"""
    (xx,yy) = ndgrid(x::AbstractVector{<:Any}, y::AbstractVector{<:Any})
todo - improve?
"""
function ndgrid(
    x::AbstractVector{<:Any},
    y::AbstractVector{<:Any},
)

    tmp = Iterators.product(x, y)
    return [p[1] for p in tmp], [p[2] for p in tmp]
end


"""
    (xx,yy,zz) = ndgrid(x::AbstractVector{<:Any}, y::..., z::...)
todo - improve?
"""
function ndgrid(
    x::AbstractVector{<:Any},
    y::AbstractVector{<:Any},
    z::AbstractVector{<:Any},
)

    tmp = Iterators.product(x, y, z)
    return [p[1] for p in tmp], [p[2] for p in tmp], [p[3] for p in tmp]
end
