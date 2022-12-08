#=
exp_xform.jl

see also MIRT/time/exp_xform.jl

2019-12-08 Connor Martin
Based on exp_xform_mex.c
Copyright 2004-9-23, Jeff Fessler, University of Michigan
=#

export exp_xform


"""
    exp_xform(x, u, v ; mode::Symbol = :matrix)

# in
* `x [N L]` possibly complex vector(s)
* `u [D N]` possibly complex vectors
* `v [D M]` possibly complex vectors
* `mode::Symbol` `:matrix` (default) | `:element` | `:row` | `:column`

# out
* `y [M L]` typically complex vector
 `y[m,l] = sum_n x[n,l] exp(-sum(u[:,n] .* v[:,m]))`

Iterates through subsets of the ML matrix designated by `:mode`
(i.e. row, column, element, or just computing the entire matrix)
This is the 'slow' 'exact' transform model for MRI.

Output type will depend on input types.
"""
function exp_xform(x::AbstractMatrix{<:Number},
        u::AbstractMatrix{<:Number},
        v::AbstractMatrix{<:Number}
        ; mode::Symbol = :matrix)

    T = promote_type(eltype(u), eltype(v), eltype(x), ComplexF32)

    N = size(u,2)
    M = size(v,2)
    out = zeros(T, M, size(x,2))

    if (mode === :matrix)
        return exp.(-(transpose(v) * u)) * x # [M D] * [D N] * [N L]

    elseif (mode === :element)
        # avoid generating intermediate steps by iterating thru N and M
        # dot the columns
        for n in 1:N
            for m in 1:M
                t = transpose(u[:,n]) * v[:,m] #calculate one spot in the array
                out[m,:] .+= (exp.(-t) .* x[n,:])
            end
        end
        return out

    elseif (mode === :row) # iterate through rows of the large matrix
        for m in 1:M
            t = transpose(v[:,m]) * u # [1 D] * [D N] -> [1 N]
            # m == 1 && @show size(t), size(exp.(-t) * x)
            # product is a 2d matrix (row), but should be a column vector:
            out[m,:] .+= (exp.(-t) * x)[1,:]
        end
        return out

    elseif (mode === :column)
        for n in 1:N
            t = transpose(v) * u[:,n] # [M D] * [D .] -> [M]
            # n == 1 && @show size(t), size(x[n,:])
            out .+= (exp.(-t) * transpose(x[n,:])) # outer product [M 1] * [1 L]
        end
        return out


    else
        throw("Invalid mode parameter $mode")
    end
end


# handle "vector" case where L=1
exp_xform(x::AbstractVector{<:Number},
        u::AbstractMatrix{<:Number},
        v::AbstractMatrix{<:Number}
        ; kwargs...) = exp_xform(reshape(x, :, 1), u, v ; kwargs...)[:,1]
