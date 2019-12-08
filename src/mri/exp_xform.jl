#=
exp_xform.jl

2019-12-08 Connor Martin
Based on exp_xform_mex.c
Copyright 2004-9-23, Jeff Fessler, University of Michigan
=#

export exp_xform

using Test: @test, @test_throws, @inferred
using BenchmarkTools: @btime


"""
    exp_xform(x, u, v ; mode::Symbol = :matrix)

in:
* `x [N L]` possibly complex vector(s)
* `u [D N]` possibly complex vectors
* `v [D M]` possibly complex vectors
* `mode::Symbol` `:matrix` (default) | `:element` | `:row` | `:column`

out:
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
        for n = 1:N
            for m = 1:M
                t = transpose(u[:,n]) * v[:,m] #calculate one spot in the array
                out[m,:] .+= (exp.(-t) .* x[n,:])
            end
        end
        return out

    elseif (mode === :row) # iterate through rows of the large matrix
        for m = 1:M
            t = transpose(v[:,m]) * u # [1 D] * [D N] -> [1 N]
            # m == 1 && @show size(t), size(exp.(-t) * x)
            # product is a 2d matrix (row), but should be a column vector:
            out[m,:] .+= (exp.(-t) * x)[1,:]
        end
        return out

    elseif (mode === :column)
        for n = 1:N
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


# test for given data type
function exp_xform_test( ; T::DataType = ComplexF32, time::Bool = false)

    modes = (:element, :row, :column)

    @info "1D tests with $T"
    N = 500
    M = 6000
    D = 3
    X = randn(T, N)
    U = randn(T, D, N)
    V = randn(T, D, M)
    y1 = @inferred exp_xform(X, U, V ; mode = :matrix)

    for mode in modes
        y2 = @inferred exp_xform(X, U, V ; mode=mode)
        @test y1 ≈ y2
    end

    @info "2D tests with $T"
    N = 10000
    M = 100
    D = 20
    L = 10
    X = randn(T, N, L)
    U = randn(T, D, N)
    V = randn(T, D, M)
    y1 = @inferred exp_xform(X, U, V ; mode = :matrix)

    for mode in modes
        y2 = @inferred exp_xform(X, U, V ; mode=mode)
        @test y1 ≈ y2
    end

    # timing tests: :matrix is fastest, with :row a close 2nd
    for mode in (:matrix, modes...)
        time && @info "Case :$mode"
        tmp = (x, u, v) -> exp_xform(x, u, v ; mode=mode)
        time && @btime $tmp($X, $U, $V)
    end

    true
end


"""
    exp_xform(:test ; time=false)
self test (with optional timing tests)
"""
function exp_xform(test::Symbol ; time::Bool = false)
    test != :test && throw("Invalid argument for exp_xform.")
	@test_throws String exp_xform(ones(2,2), ones(2,2), ones(2,2) ; mode=:bad)

#    for T in (ComplexF32, ComplexF64)
        @test exp_xform_test( ; time=time) # T=T
#    end
    true
end

#exp_xform(:test ; time=true)
