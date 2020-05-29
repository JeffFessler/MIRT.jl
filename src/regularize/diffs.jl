#=
diffs.jl
finite differences
2019-03-06 Jeff Fessler, University of Michigan
=#

export diff_map

using LinearMapsAA: LinearMapAA
using Test: @test, @test_throws
using BenchmarkTools: @btime


"""
`d = diff2d_forw(X)`

2D finite differences along both dimensions, for anisotropic TV regularization.
Performs the same operations as
d = ``[(I_N \\otimes D_M); (D_N \\otimes I_M)] X[:]``
where ``D_N`` denotes the ``N-1 \\times N`` 1D finite difference matrix
and ``\\otimes`` denotes the Kronecker product,
but does it efficiently
without using `spdiagm` (or any `SparseArrays` function).

in
- `X`		`M x N` array (typically a 2D image).
It cannot be a Vector!  (But it can be a `Mx1` or `1xN` 2D array.)

out
- `d`		vector of length `N*(M-1) + (N-1)*M`
"""
function diff2d_forw(x::AbstractMatrix{<:Number})
	return [diff(x,dims=1)[:]; diff(x,dims=2)[:]]
end


"""
`d = diffnd_forw(X)`

N-D finite differences along all dimensions, for anisotropic TV regularization.
Performs the same operations as
d = ``[(I_{N_d} \\otimes \\cdots \\otimes D_{N_1}); \\dots; (D_{N_d} \\otimes \\cdots \\otimes I_{N_1})] X[:]``
where ``D_N`` denotes the ``N-1 \\times N`` 1D finite difference matrix
and ``\\otimes`` denotes the Kronecker product,
but does it efficiently
without using `spdiagm` (or any `SparseArrays` function).

in
- `X`		`N_1 x ... x N_d` array (typically an N-D image).

out
- `d`		vector of length `N_d*...*(N_1-1) + ... + (N_d-1)*...*N_1`
"""
function diffnd_forw(x::AbstractArray{<:Number,D}) where {D}
    return reduce(vcat, vec(diff(x, dims = d)) for d = 1:D)
end


"""
`z =  diff2d_adj(d, M, N; out2d=false)`

Adjoint of 2D finite differences along both dimensions.
Performs the same operations as
``z = [(I_N \\otimes D_M); (D_N \\otimes I_M)]' * d``
where D_N denotes the N-1 x N 1D finite difference matrix
and \\otimes denotes the Kronecker product,
but does it efficiently without using spdiagm (or any SparseArrays function).

in
- `d`		vector of length `N*(M-1) + (N-1)*M`
- `M,N`		desired output size

option
- `out2d`	if true then return `M x N` array, else `M*N` vector

out
- `z`		`M*N` vector or `M x N` array (typically a 2D image)

"""
function diff2d_adj(d::AbstractVector{<:Number}, M::Int, N::Int ; out2d=false)

	length(d) != N*(M-1) + (N-1)*M && throw("length(d)")

	dx = d[1:(N*(M-1))]
	dx = reshape(dx, M-1, N)
	dy = d[1+(N*(M-1)):end]
	dy = reshape(dy, M, N-1)
	zx = [-transpose(dx[1,:]);
		(@views dx[1:end-1,:] - dx[2:end,:]);
		transpose(dx[end,:])]
	zy = [-dy[:,1] (@views dy[:,1:end-1] - dy[:,2:end]) dy[:,end]]
	z = zx + zy # M by N

	return out2d ? z : z[:] # array or vector?
end


"""
`z =  diffnd_adj(d, N...; out2d=false)`

Adjoint of N-D finite differences along both dimensions.
Performs the same operations as
``z = [(I_{N_d} \\otimes \\cdots \\otimes D_{N_1}); \\dots; (D_{N_d} \\otimes \\cdots \\otimes I_{N_1})]' * d``
where D_N denotes the N-1 x N 1D finite difference matrix
and \\otimes denotes the Kronecker product,
but does it efficiently without using spdiagm (or any SparseArrays function).

in
- `d`		vector of length `N_d*...*(N_1-1) + ... + (N_d-1)*...*N_1`
- `N...`	desired output size

option
- `outnd`	if true then return `N_1 x ... x N_d` array, else `prod(N)` vector

out
- `z`		`prod(N)` vector or `N_1 x ... x N_d` array (typically an N-D image)

"""
function diffnd_adj(d::AbstractVector{<:Number}, N::Int... ; outnd=false)

    # Note that N must be strictly greater than 1 for each dimension,
    # or N must be 1 for all dimensions
    # (This is true of diff2d_adj as well)

    ndims = length(N)
    length(d) != sum(*(N[1:i-1]..., N[i] - 1, N[i+1:end]...) for i = 1:ndims) &&
        throw("length(d)")

    z = zeros(eltype(d), N...)
    for i = 1:ndims
        if i == 1
            di = @view(d[1:*(N[i] - 1, N[i+1:end]...)])
        else
            start = 1 + sum(*(N[1:i-n-1]..., N[i-n] - 1, N[i-n+1:end]...) for n = 1:i-1)
            len = *(N[1:i-1]..., N[i] - 1, N[i+1:end]...)
            di = @view(d[start:start+len-1])
        end
        di = reshape(di, N[1:i-1]..., N[i] - 1, N[i+1:end]...)
        slice1 = selectdim(z, i, 1)
        slice1 .-= selectdim(di, i, 1)
        slicen = selectdim(z, i, 2:N[i]-1)
        slicen .+= selectdim(di, i, 1:N[i]-2) - selectdim(di, i, 2:N[i]-1)
        sliceN = selectdim(z, i, N[i])
        sliceN .+= selectdim(di, i, N[i] - 1)
    end

    return outnd ? z : vec(z)
end


"""
`T = diff2d_map(M::Int, N::Int)`
"""
function diff2d_map(M::Int, N::Int)
	return LinearMapAA(
        x -> diff2d_forw(reshape(x,M,N)),
        d -> diff2d_adj(d, M, N),
        (N*(M-1)+M*(N-1), N*M), (name="diff2_map",))
end


"""
`T = diffnd_map(N::Int...)`
"""
function diffnd_map(N::Int...)
    return LinearMapAA(
    x -> diffnd_forw(reshape(x,N...)),
    d -> diffnd_adj(d, N...),
    (sum(*(N[1:i-1]..., N[i] - 1, N[i+1:end]...) for i = 1:length(N)), prod(N)),
    (name="diffn_map",))
end


"""
`T = diff_map(M::Int, N::Int)`

in
- `M,N` image size

out
- `T` a `LinearMapAA` object for regularizing via `T*x`
"""
function diff_map(M::Int, N::Int)
	return diff2d_map(M,N)
end


"""
`diff_map(:test)`
self test
"""
function diff_map(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	M,N = 4,5
	T = diff_map(M,N)
	@test Matrix(T)' == Matrix(T') # adjoint test
	@test T.name == "diff2_map"
	true
end


"""
`diffnd_map(:test)`
self test
"""
function diffnd_map(test::Symbol)
    test != :test && throw(ArgumentError("test $test"))
    for N in [(2,), (10,), (2,3), (10,11), (1,1,1), (2,3,4), (4,4,4,4)]
        T = diffnd_map(N...)
        @test Matrix(T)' == Matrix(T')
        @test T.name == "diffn_map"
    end
    # adjoint doesn't work if any of the dimensions has size 1
    # (unless all are size 1)
    N = (1,2)
    T = diffnd_map(N...)
    @test_throws BoundsError Matrix(T)' == Matrix(T')
    # compare to diff_map
    N = (100, 101)
    T2d = diff_map(N...)
    Tnd = diffnd_map(N...)
    x = randn(prod(N))
    d2d = @btime $T2d * $x # 38 μs
    dnd = @btime $Tnd * $x # 27 μs
    @test d2d == dnd
    d = d2d
    y2d = @btime $T2d' * $d # 65 μs
    ynd = @btime $Tnd' * $d # 60 μs
    @test y2d == ynd
    true
end
