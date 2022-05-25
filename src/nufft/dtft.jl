#=
dtft.jl
discrete-time Fourier transform
The main purpose of this code is for testing the NUFFT code.
2019-06-05, Jeff Fessler, University of Michigan
Note: the multi-dimensional version uses collect(Tuple(idx))
that might be inefficient, but for efficiency one should just use NUFFT.

see also MIRT/time/dtft.jl
=#

export dtft_init, dtft, dtft_adj

#using LinearAlgebra: norm
using Distributed: @sync, @distributed, pmap
using SharedArrays: SharedArray, sdata
using LinearMapsAA: LinearMapAA, LinearMapAO


"""
    X = dtft(w, x ; n_shift=?)
1D DTFT

``X[m] = \\sum_{n=0}^{N-1} x[n] \\exp(-i w[m] (n - n_{shift})), m=1,…,M``

in
- `w::AbstractVector{<:Real}`	`[M]` frequency locations ("units" radians/sample)
- `x::AbstractVector{<:Number}`	`[N]` 1D signal

option
- `n_shift::Real` often is N/2; default 0

out
- `X::AbstractVector{ComplexF64}` `[M]` DTFT
"""
function dtft(w::AbstractVector{<:Real}, x::AbstractVector{<:Number}
		; n_shift::Real = 0)
	return dtft_loop_n(w, x ; n_shift=n_shift)
end


"""
    X = dtft(w, x ; n_shift=?)
multi-dimensional DTFT (DSFT)

``X[m] = \\sum_{n=0}^{N-1} x[n] \\exp(-i w[m,:] (n - n_{shift})), m=1,…,M``
where here `n` is a `CartesianIndex`

in
- `w::AbstractMatrix{<:Real}`	`[M,D]` frequency locations ("units" radians/sample)
- `x::AbstractArray{<:Number}`	`[(Nd)]` D-dimensional signal

option
- `n_shift::AbstractVector{<:Real}` often is N/2; default zeros(D)

out
- `X::AbstractVector{ComplexF64}` `[M]` DTFT
"""
function dtft(w::AbstractMatrix{<:Real}, x::AbstractMatrix{<:Number}
	; n_shift::AbstractVector{<:Real} = zeros(Int, ndims(x)),
)
	return dtft_loop_n(w, x ; n_shift=n_shift)
end



"""
    d = dtft_init(w, N ; n_shift=?)
for 1D DTFT

in
- `w::AbstractVector{<:Real}`	`[M]` frequency locations ("units" radians/sample)
- `N::Int` size of signal `x`

option
- `n_shift::Real` often is N/2; default 0
- `T::DataType` default `ComplexF64` for testing NUFFT accuracy

out
- `d::NamedTuple` with fields
	`dtft = x -> dtft(x), adjoint = y -> dtft_adj(y), A=LinearMapAO`
"""
function dtft_init(
	w::AbstractVector{<:Real},
	N::Int ;
	n_shift::Real = 0,
	T::DataType = ComplexF64,
)

	M = length(w)
	forw = x -> dtft(w, x ; n_shift=n_shift)
	back = y -> dtft_adj(w, y, N ; n_shift=n_shift)
	A = LinearMapAA(forw, back, (M, N),
		(name="dtft1", N=(N,)) ; T = ComplexF64,
		operator = true,
	)
	return (dtft=forw, adjoint=back, A=A)
end


"""
    d = dtft_init(w, N ; n_shift=?)
for multi-dimensional DTFT (DSFT)

in
- `w::AbstractMatrix{<:Real}`	`[M,D]` frequency locations ("units" radians/sample)
- `N::Dims` `[D]` dimensions of signal `x`

option
- `n_shift::AbstractVector{<:Real}` often is N/2; default zeros(D)
- `T::DataType` default `ComplexF64` for testing NUFFT accuracy

out
- `d::NamedTuple` with fields
	`dtft = x -> dtft(x), adjoint = y -> dtft_adj(y), A=LinearMapAO`
"""
function dtft_init(
	w::AbstractMatrix{<:Real},
	N::Dims ;
	n_shift::AbstractVector{<:Real} = zeros(Int, length(N)),
)

	D = size(w,2)
	length(N) != D && throw(DimensionMismatch("length(N) vs D=$D"))
	length(n_shift) != D && throw(DimensionMismatch("length(n_shift) vs D=$D"))
	M = size(w,1)
	forw = x -> dtft(w, x ; n_shift=n_shift)
	back = y -> dtft_adj(w, y, N ; n_shift=n_shift)
	A = LinearMapAA(forw, back, (M, prod(N)),
		(name="dtft$(length(N))", N=N) ; T = ComplexF64,
		operator = true, idim = N,
	)
	return (dtft=forw, adjoint=back, A=A)
end


"""
    x = dtft_adj(w, X, N ; n_shift=?)
adjoint for 1D DTFT

``x[n] = \\sum_{m=1}^M X[m] \\exp(i w[m] (n - n_{shift})), n=0,…,N-1``

This is the *adjoint* (transpose) of `dtft`, not an *inverse* DTFT.

in
- `w::AbstractVector{<:Real}`	`[M]` frequency locations ("units" radians/sample)
- `X::AbstractVector{ComplexF64}`	`[M]` spectrum values
- `N::Int` size of signal `x`

option
- `n_shift::Real` often is N/2; default 0
- `T::DataType` output data type; default `ComplexF64`

out
- `x::AbstractVector{<:Number}`	signal [N]
"""
function dtft_adj(w::AbstractVector{<:Real}, X::AbstractVector{<:Number},
	N::Int ; n_shift::Real = 0,
	T::DataType = ComplexF64,
)
	M = length(w)
	out = similar(X, T, N)
	nshift1 = n_shift + 1
	X = T.(X)
	for n=1:N # loop over output signal samples
		out[n] = sum(cis.((n - nshift1) * w) .* X)
	end
	return out
end


"""
    x = dtft_adj(w, X, N ; n_shift=?)
adjoint for multi-dimensional DTFT (DSFT)

``x[n] = \\sum_{m=1}^M X[m] \\exp(i w[m,:] (n - n_{shift})), n=0,…,N-1``
where here `n` is a `CartesianIndex`

in
- `X::AbstractVector{ComplexF64}` `[M]` DTFT
- `w::AbstractMatrix{<:Real}` `[M,D]` frequency locations ("units" radians/sample)
- `N::Dims`	`[D]` dimensions of signal `x`

option
- `n_shift::AbstractVector{<:Real}` often is `N/2`; default `zeros(D)`
- `T::DataType` default `(eltype(w) == Float64) ? ComplexF64 : ComplexF32`

out
- `x::AbstractArray{<:Number}` `[(N)]` `D`-dimensional signal
"""
function dtft_adj(
	w::AbstractMatrix{<:Real},
	X::AbstractVector{<:Number},
	N::Dims ;
	n_shift::AbstractVector{<:Real} = zeros(Int, ndims(x)),
	T::DataType = (eltype(w) == Float64) ? ComplexF64 : ComplexF32,
)

	M, D = size(w)
	nshift1 = n_shift .+ 1

	out = similar(X, T, N)

	for n=1:prod(N)
		idx = CartesianIndices(N)[n]
		tmp = cis.(w * (collect(Tuple(idx)) - nshift1))
		tmp = transpose(cis.(w * (collect(Tuple(idx)) - nshift1))) * X
		out[idx] = cis.(-w * (collect(Tuple(idx)) - nshift1))' * X
	end
	return out
end


# 1D loop over signal values
function dtft_loop_n(w::AbstractVector{<:Real}, x::AbstractVector{<:Number}
	; n_shift::Real = 0,
)
	N = length(x)
	out = fill(ComplexF64(0), size(w))
	nshift1 = n_shift + 1
	for n=1:N # tried @simd but not faster
		out .+= x[n] * cis.(-(n-nshift1) * w)
	end
	return out
end


# dD loop over signal values
# * `w [M,D]`
# * `x [(Nd)]`
function dtft_loop_n(
	w::AbstractMatrix{<:Real},
	x::AbstractMatrix{<:Number} ;
	n_shift::AbstractVector{<:Real} = zeros(Int, size(w,2)),
)
	M,D = size(w)
	N = size(x)
	length(N) != D && throw(DimensionMismatch("length(N) vs D=$D"))
	length(n_shift) != D && throw(DimensionMismatch("length(n_shift) vs D=$D"))
	ndims(x) != D && throw(DimensionMismatch("ndims(x) vs D=$D"))
	nshift1 = n_shift .+ 1
	T = (eltype(w) == Float64) ? ComplexF64 : ComplexF32
	out = fill(T(0), M)
	for n=1:prod(N)
		idx = CartesianIndices(N)[n]
		out .+= x[idx] * cis.(-w * (collect(Tuple(idx)) - nshift1))
	end
	return out
end


# 1D loop over frequencies
function dtft_loop_m(
	w::AbstractVector{<:Real},
	x::AbstractVector{<:Number} ;
	n_shift::Real = 0,
    T::DataType = ComplexF64,
)
	M = length(w)
	N = length(x)
	out = similar(x, T, M)
	nn = -((0:(N-1)) .- n_shift)
	x = T.(x)
	for m=1:M
		out[m] = sum(cis.(nn * w[m]) .* x)
	end
	return out
end


# 1D matrix * vector
function dtft_matvec(
	w::AbstractVector{<:Real},
	x::AbstractVector{<:Number} ;
	n_shift::Real = 0,
)
	N = length(x)
	return cis.(w * (-((0:(N-1)) .- n_shift))') * x
end


# 1D single frequency
function dtft_one_w(
	w::Real,
	x::AbstractVector{<:Number} ;
	n_shift::Real = 0,
)
	N = length(x)
	return sum(cis.(-((0:(N-1)) .- n_shift) * w) .* x)
end


# 1D pmap
function dtft_pmap_m(
	w::AbstractVector{<:Real},
	x::AbstractVector{<:Number}
	; n_shift::Real = 0,
)
	tmp = w -> dtft_one_w(w, x ; n_shift=n_shift)
	return pmap(tmp, w)
end


# 1D distributed loop over frequencies
function dtft_dist_m(
	w::AbstractVector{<:Real},
	x::AbstractVector{<:Number} ;
	n_shift::Real = 0,
)
	M = length(w)
	N = length(x)
	out = SharedArray{ComplexF64}(M)
	nn = -((0:(N-1)) .- n_shift)
	x = ComplexF64.(x)
	@sync @distributed for m=1:M
		out[m] = sum(cis.(nn * w[m]) .* x)
	end
	return sdata(out)
end
