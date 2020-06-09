#=
diffs.jl
finite differences
2019-03-06 Jeff Fessler, University of Michigan
2020-06 N-D version by Steven Whitaker
=#

export diff_map

using LinearMapsAA: LinearMapAA
using Test: @test, @test_throws


# size after apply `diff` along `dim`
diff_size(N::Dims, dim::Int) = (N[1:dim-1]..., N[dim] - 1, N[dim+1:end]...)

# corresponding length
diff_length(N::Dims, dim::Int) = prod(diff_size(N, dim))

# check argument validity
function diff_check(N::Dims, dims::AbstractVector{Int})
    !all(1 .<= dims .<= length(N)) && throw(ArgumentError("dims range $dims"))
    (length(dims) > 1) && (unique(dims) != dims) &&
        throw(ArgumentError("non-unique dims $dims"))
    any(N[dims] .== 1) && throw(ArgumentError("invalid size $N for dims $dims"))
end


"""
    d = diff_forw(X ; dims::AbstractVector{Int} = 1:ndims(X))

Finite differences along one or more dimensions of an array,
e.g., for anisotropic TV regularization.

By default performs the same operations as
``d = [(I_{N_d} \\otimes \\cdots \\otimes D_{N_1}); \\dots; (D_{N_d} \\otimes \\cdots \\otimes I_{N_1})] X[:]``
where ``D_N`` denotes the `N-1 × N` 1D finite difference matrix
and ``\\otimes`` denotes the Kronecker product, but does it efficiently
without using `spdiagm` (or any `SparseArrays` function).

Input dimension `N` must exceed 1 for each dimension specified by `dims`.

in
- `X` `N_1 × ... × N_d` array (typically an N-D image).

option
- `dims` dimensions along which to perform finite differences; default `1:ndims(X)`
must have unique elements and be a subset of 1:ndims(X)

out
- `d` vector of default length `N_d*...*(N_1-1) + ... + (N_d-1)*...*N_1`
"""
function diff_forw(x::AbstractArray{<:Number,D}
    ; dims::AbstractVector{Int} = 1:D,
) where {D}
    diff_check(size(x), dims)
    return reduce(vcat, vec(diff(x, dims = d)) for d in dims)
end


#=
I wanted to give flexibility for Ints and Tuples
but I could not get this to work...

"""
    d = diff_forw(X ; dims::AbstractVector{Int} = 1:ndims(X))
"""
diff_forw(x::AbstractArray{<:Number} ; dims::AbstractVector{Int} = 1:ndims(X)) =
    diff_forw(x, dims)

"""
    d = diff_forw(X ; dims::Dims = (1,))
"""
diff_forw(x::AbstractArray{<:Number} ; dims::Dims = (1,)) =
    diff_forw(x, collect(dims))

"""
    d = diff_forw(X ; dims::Int = 1)
"""
diff_forw(x::AbstractArray{<:Number} ; dims::Int = 1) =
    diff_forw(x, [dims])
=#


"""
    Z = diff_adj(dx, N... ; dims=1:length(N))

Adjoint of finite differences of arrays along one or more dimensions.
By default performs the same operations as
``vec(Z) = [(I_{N_d} \\otimes \\cdots \\otimes D_{N_1}); \\dots; (D_{N_d} \\otimes \\cdots \\otimes I_{N_1})]' * d``
where D_N denotes the `N-1 × N` 1D finite difference matrix
and `\\otimes` denotes the Kronecker product,
but does it efficiently without using `spdiagm` (or any `SparseArrays` function).

in
- `dx` vector of typical length `N_d*...*(N_1-1) + ... + (N_d-1)*...*N_1`
- `N...` desired output size

option
- `dims` dimensions along which to perform adjoint finite differences

out
- `Z` `N_1 × ... × N_d` array by default

"""
function diff_adj(d::AbstractVector{<:Number}, N::Vararg{Int,D}
    ; dims::AbstractVector{Int} = 1:D,
) where {D}

# todo: N::Dims instead of Vararg?

    length(d) != sum(diff_length(N,dim) for dim in dims) && throw("length(d)")

    z = zeros(eltype(d), N...)
    for (i, dim) in enumerate(dims)
        if i == 1
            di = @view(d[1:diff_length(N,dim)])
        else
            start = 1 + sum(diff_length(N,n) for n in dims[1:i-1])
            len = diff_length(N,dim)
            di = @view(d[start:start+len-1])
        end
        di = reshape(di, diff_size(N,dim))
        slice1 = selectdim(z, dim, 1)
        slice1 .-= selectdim(di, dim, 1)
        slicen = selectdim(z, dim, 2:N[dim]-1)
        slicen .+= selectdim(di, dim, 1:N[dim]-2) - selectdim(di, dim, 2:N[dim]-1)
        sliceN = selectdim(z, dim, N[dim])
        sliceN .+= selectdim(di, dim, N[dim] - 1)
    end

    return z
end

# backward compatibility, undocumented because deprecated
function diff2d_forw(x::AbstractMatrix{<:Number})
    isinteractive() && @warn("diff2d_forw is deprecated; use diff_forw")
    return diff_forw(x)
end

function diff2d_adj(d::AbstractVector{<:Number}, M::Int, N::Int ; out2d=false)
    isinteractive() && @warn("diff2d_adj is deprecated; use diff_adj")
    tmp = diff_adj(d, M, N)
    return out2d ? tmp : tmp[:]
end

function diff2d_map(M::Int, N::Int)
    isinteractive() && @warn("diff2d_map is deprecated; use diff_map")
    return diff_map(M, N)
end


"""
    T = diff_map(N::Int... ; dims=1:length(N))

in
- `N...` image size

out
- `T` `LinearMapAA` object for regularizing via `T*x`
"""
function diff_map(N::Vararg{Int,D} ; dims::AbstractVector{Int} = 1:D) where {D}
    diff_check((N...,), dims)
    return LinearMapAA(
        x -> diff_forw(reshape(x, N...), dims=dims),
        d -> diff_adj(d, N..., dims=dims)[:],
        (sum(diff_length(N,dim) for dim in dims), prod(N)),
        (name="diff_map", dims=dims),
    )
end


"""
`diff_map(:test)`
self test
"""
function diff_map(test::Symbol)
    test != :test && throw(ArgumentError("test $test"))
    @test_throws ArgumentError diff_forw(ones(3), dims=[2])
    @test_throws ArgumentError diff_forw(ones(3), dims=[1,2])
    @test_throws ArgumentError diff_forw(ones(1,2,1), dims=[1,2])
    @test diff_forw(ones(2,4,6) ; dims=[2]) == zeros(2*3*6)

    for N in ([2,], [10,], [2,3], [10,11], [2,3,4], [4,4,4,4]) # [1,1,1],
        T = diff_map(N...)
        @test Matrix(T)' == Matrix(T')
        @test T.name == "diff_map"
        T = diff_map(N..., dims=[1])
        @test Matrix(T)' == Matrix(T')
        if length(N) >= 2
            for dims in [[2], [1,2]]
                T = diff_map(N..., dims=dims)
                @test Matrix(T)' == Matrix(T')
            end
        end
        if length(N) >= 3
            for dims in [[3], [1,3], [2,3], [1,2,3]]
                T = diff_map(N..., dims=dims)
                @test Matrix(T)' == Matrix(T')
            end
        end
    end

    if true # test old 2D versions
        N = (3,5)
        x = rand(N...)
        @test diff2d_forw(x) == diff_forw(x)
        d = rand(sum(diff_length(N,dim) for dim=1:2))
        @test diff2d_adj(d, N... ; out2d=false) == diff_adj(d, N...)[:]
        @test diff2d_adj(d, N... ; out2d=true) == diff_adj(d, N...)
        @test Matrix(diff2d_map(N...)) == Matrix(diff_map(N...))
    end

    # adjoint doesn't work if any of the dimensions specified by dims has size 1
    # (unless all are size 1)
    N = (1,2)
    @test_throws ArgumentError T = diff_map(N...)
    @test_throws ArgumentError T = diff_map(N..., dims = [1])
    T = diff_map(N..., dims = [2])
    @test Matrix(T)' == Matrix(T')
    true
end

# diff_map(:test)
