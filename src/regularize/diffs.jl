# diffs.jl
# finite differences
# 2019-03-96 Jeff Fessler, University of Michigan

using LinearMaps

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
* `X`		`M x N` array (typically a 2D image).
It cannot be a Vector!  (But it can be a `Mx1` or `1xN` 2D array.)

out
* `d`		vector of length `N*(M-1) + (N-1)*M`
"""
function diff2d_forw(x::AbstractMatrix{<:Number})

return [diff(x,dims=1)[:]; diff(x,dims=2)[:]]
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
* `d`		vector of length `N*(M-1) + (N-1)*M`
* `M,N`		desired output size

option
* `out2d`	if true then return `M x N` array, else `M*N` vector

out
* `z`		`M*N` vector or `M x N` array (typically a 2D image)

"""
function diff2d_adj(d::AbstractVector{<:Number}, M, N; out2d=false)

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
`T = diff2d_map(M,N)`
"""
function diff2d_map(M,N)
	return LinearMap(
        x -> diff2d_forw(reshape(x,M,N)),
        d -> diff2d_adj(d, M, N),
        N*(M-1)+M*(N-1), N*M)
end


"""
`T = diff_map(M,N)`

in
* `M,N` image size

out
* `T` a `LinearMap` object for regularizing via T*x
"""
function diff_map(M,N)
	return diff2d_map(M,N)
end


"""
`diff_map(:test)` self test
"""
function diff_map(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	M,N = 4,5
	T = diff_map(M,N)
	@test Matrix(T)' == Matrix(T') # adjoint test
	true
end
