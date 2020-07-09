#=
reverser.jl
based on
https://stackoverflow.com/questions/27411401/julia-reverse-n-dimensional-arrays
=#

export reverser


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
