#=
dtft.jl
discrete-time Fourier transform
The main purpose of this code is for testing the NUFFT code.
2019-06-05, Jeff Fessler, University of Michigan
Note: the multi-dimensional version uses collect(Tuple(idx))
that might be inefficient, but for efficiency one should just use NUFFT.
=#

export dtft_init, dtft, dtft_adj


using FFTW: fft
using Random: seed!
#using LinearAlgebra: norm
#using BenchmarkTools: @btime
using Distributed: @sync, @distributed, pmap
using SharedArrays: SharedArray, sdata
using Test: @test
using LinearMaps: LinearMap


"""
`X = dtft(w, x; n_shift=?)` 1D DTFT

``X[m] = sum_{n=0}^{N-1} x[n] exp(-i w[m] (n - n_shift)), m=1,…,M``

in
- `w::AbstractVector{<:Real}`	`[M]` frequency locations ("units" radians/sample)
- `x::AbstractVector{<:Number}`	`[N]` 1D signal

option
- `n_shift::Real` often is N/2; default 0

out
- `X::AbstractVector{ComplexF64}`		`[M]` DTFT
"""
function dtft(w::AbstractVector{<:Real}, x::AbstractVector{<:Number}
		; n_shift::Real = 0)
	return dtft_loop_n(w, x; n_shift=n_shift)
end


"""
`X = dtft(w, x; n_shift=?)` multi-dimensional DTFT (DSFT)

``X[m] = sum_{n=0}^{N-1} x[n] exp(-i w[m,:] (n - n_shift)), m=1,…,M``
where here `n` is a `CartesianIndex`

in
- `w::AbstractMatrix{<:Real}`	`[M,D]` frequency locations ("units" radians/sample)
- `x::AbstractArray{<:Number}`	`[(Nd)]` D-dimensional signal

option
- `n_shift::AbstractVector{<:Real}` often is N/2; default zeros(D)

out
- `X::AbstractVector{ComplexF64}`		`[M]` DTFT
"""
function dtft(w::AbstractMatrix{<:Real}, x::AbstractMatrix{<:Number};
	n_shift::AbstractVector{<:Real} = zeros(Int, ndims(x)))
	return dtft_loop_n(w, x; n_shift=n_shift)
end



"""
`d = dtft_init(w, N; n_shift=?)` for 1D DTFT

in
- `w::AbstractVector{<:Real}`	`[M]` frequency locations ("units" radians/sample)
- `N::Int` size of signal `x`

option
- `n_shift::Real` often is N/2; default 0

out
- `d::NamedTuple` with fields
	`dtft = x -> dtft(x), adjoint = y -> dtft_adj(y), A=LinearMap`
"""
function dtft_init(w::AbstractVector{<:Real}, N::Int; n_shift::Real = 0)
	M = length(w)
	forw = x -> dtft(w, x; n_shift=n_shift)
	back = y -> dtft_adj(w, y, N; n_shift=n_shift)
	A = LinearMap{ComplexF64}(x -> forw(x), y -> back(y), M, N)
	return (dtft=forw, adjoint=back, A=A)
end


"""
`d = dtft_init(w, N; n_shift=?)` for multi-dimensional DTFT (DSFT)

in
- `w::AbstractMatrix{<:Real}`	`[M,D]` frequency locations ("units" radians/sample)
- `N::Dims` `[D]` dimensions of signal `x`

option
- `n_shift::AbstractVector{<:Real}` often is N/2; default zeros(D)

out
- `d::NamedTuple` with fields
	`dtft = x -> dtft(x), adjoint = y -> dtft_adj(y), A=LinearMap`
"""
function dtft_init(w::AbstractMatrix{<:Real}, N::Dims;
		n_shift::AbstractVector{<:Real} = zeros(Int, length(N)))
	D = size(w,2)
	length(N) != D && throw(DimensionMismatch("length(N) vs D=$D"))
	length(n_shift) != D && throw(DimensionMismatch("length(n_shift) vs D=$D"))
	M = size(w,1)
	forw = x -> dtft(w, reshape(x, N); n_shift=n_shift)
	back = y -> dtft_adj(w, y, N; n_shift=n_shift)
	A = LinearMap{ComplexF64}(x -> forw(x), y -> back(y)[:], M, prod(N))
	return (dtft=forw, adjoint=back, A=A)
end


"""
`x = dtft_adj(w, X, N; n_shift=?)` adjoint for 1D DTFT

``x[n] = sum_{m=1}^M X[m] exp(i w[m] (n - n_shift)), n=0,…,N-1``

This is the *adjoint* (transpose) of `dtft`, not an *inverse* DTFT.

in
- `w::AbstractVector{<:Real}`	`[M]` frequency locations ("units" radians/sample)
- `X::AbstractVector{ComplexF64}`	`[M]` spectrum values
- `N::Int` size of signal `x`

option
- `n_shift::Real` often is N/2; default 0

out
- `x::AbstractVector{<:Number}`	signal [N]
"""
function dtft_adj(w::AbstractVector{<:Real}, X::AbstractVector{<:Number},
		N::Int; n_shift::Real = 0)
	M = length(w)
	out = Array{ComplexF64}(undef, N)
	nshift1 = n_shift + 1
	X = ComplexF64.(X)
	for n=1:N # loop over output signal samples
		out[n] = sum(cis.((n - nshift1) * w) .* X)
	end
	return out
end


"""
`x = dtft_adj(w, X, N; n_shift=?)` adjoint for multi-dimensional DTFT (DSFT)

``x[n] = sum_{m=1}^M X[m] exp(i w[m,:] (n - n_shift)), n=0,…,N-1``
where here `n` is a `CartesianIndex`

in
- `X::AbstractVector{ComplexF64}`		`[M]` DTFT
- `w::AbstractMatrix{<:Real}`	`[M,D]` frequency locations ("units" radians/sample)
- `N::Dims`	`[D]` dimensions of signal `x`

option
- `n_shift::AbstractVector{<:Real}` often is N/2; default zeros(D)

out
- `x::AbstractArray{<:Number}`	`[(N)]` D-dimensional signal
"""
function dtft_adj(w::AbstractMatrix{<:Real}, X::AbstractVector{<:Number},
		N::Dims ;
		n_shift::AbstractVector{<:Real} = zeros(Int, ndims(x)))

	M,D = size(w)
	nshift1 = n_shift .+ 1

	T = (eltype(w) == Float64) ? ComplexF64 : ComplexF32
	out = Array{T}(undef, N)

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
		; n_shift::Real = 0)
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
function dtft_loop_n(w::AbstractMatrix{<:Real}, x::AbstractMatrix{<:Number}
		; n_shift::AbstractVector{<:Real} = zeros(Int, size(w,2)))
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
function dtft_loop_m(w::AbstractVector{<:Real}, x::AbstractVector{<:Number}
		; n_shift::Real = 0)
	M = length(w)
	N = length(x)
	out = Array{ComplexF64}(undef, M)
	nn = -((0:(N-1)) .- n_shift)
	x = ComplexF64.(x)
	for m=1:M
		out[m] = sum(cis.(nn * w[m]) .* x)
	end
	return out
end


# 1D matrix * vector
function dtft_matvec(w::AbstractVector{<:Real}, x::AbstractVector{<:Number}
		; n_shift::Real = 0)
	N = length(x)
	return cis.(w * (-((0:(N-1)) .- n_shift))') * x
end


# 1D single frequency
function dtft_one_w(w::Real, x::AbstractVector{<:Number}
		; n_shift::Real = 0)
	N = length(x)
	return sum(cis.(-((0:(N-1)) .- n_shift) * w) .* x)
end


# 1D pmap
function dtft_pmap_m(w::AbstractVector{<:Real}, x::AbstractVector{<:Number}
		; n_shift::Real = 0)
	tmp = w -> dtft_one_w(w, x; n_shift=n_shift)
	return pmap(tmp, w)
end


# 1D distributed loop over frequencies
function dtft_dist_m(w::AbstractVector{<:Real}, x::AbstractVector{<:Number}
		; n_shift::Real = 0)
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


# 1D small test data to verify correctness
function dtft_test1( ; N::Int=2^10)
	seed!(0)
	x = randn(ComplexF64, N)
	w = (2*pi) * (0:N-1) / N

	o0 = fft(x)

	o1 = dtft_loop_n(w, x)
	o2 = dtft_loop_m(w, x)
	o3 = dtft_dist_m(w, x)
	o4 = dtft_pmap_m(w, x)
	o5 = dtft_matvec(w, x)
#	@show norm(o1 - o0) / norm(o0)

	@test isapprox(o1, o0)
	@test isapprox(o2, o0)
	@test isapprox(o3, o0)
	@test isapprox(o4, o0)
	@test isapprox(o5, o0)

	true
end


# 2D small test
function dtft_test2( ; N::Dims = (2^3,2^2))
	seed!(0)
	x = randn(ComplexF64, N)
	w1 = (2*pi) * (0:N[1]-1) / N[1]
	w2 = (2*pi) * (0:N[2]-1) / N[2]
	w1 = repeat(w1, 1, N[2])
	w2 = repeat(w2', N[1], 1)
	w = [w1[:] w2[:]]
	d = dtft_init(w, N)

	o0 = fft(x)
#	o1 = dtft_loop_n(w, x)
	o1 = d.dtft(x)
	o1 = reshape(o1, N)
	@test isapprox(o1, o0)
	true
end


# 1D verify consistency
function dtft_test1c( ; N::Int=2^10, M::Int=2^11, n_shift::Real=7)
	seed!(0)
	x = randn(ComplexF64, N)
	w = randn(M)

	o1 = dtft_loop_n(w, x; n_shift=n_shift)
	o2 = dtft_loop_m(w, x; n_shift=n_shift)
	o3 = dtft_dist_m(w, x; n_shift=n_shift)
	o4 = dtft_pmap_m(w, x; n_shift=n_shift)
	o5 = dtft_matvec(w, x; n_shift=n_shift)

	@test isapprox(o2, o1)
	@test isapprox(o3, o1)
	@test isapprox(o4, o1)
	@test isapprox(o5, o1)

	d = dtft_init(w, N; n_shift=n_shift)
	o6 = d.dtft(x)
	@test isapprox(o6, o1)

	o7 = d.A * x
	@test isapprox(o7, o1)

	b1 = dtft_adj(w, o1, N; n_shift=n_shift)
	b2 = d.adjoint(o1)
	@test isapprox(b2, b1)
	b3 = d.A' * o1
	@test isapprox(b3, b1)

#=
	# time DTFT
	@btime fft($x)
	@btime dtft_loop_n($w, $x)
	@btime dtft_loop_m($w, $x)
	@btime dtft_dist_m($w, $x)
	@btime dtft_pmap_m($w, $x)
	@btime dtft_matvec($w, $x)

	# jf28 results:
	# 12.351 μs (51 allocations: 19.08 KiB)
	# 38.608 ms (6146 allocations: 80.34 MiB)
	# 43.119 ms (2051 allocations: 32.30 MiB)
	# 44.179 ms (11965 allocations: 32.64 MiB)
	# 51.865 ms (20518 allocations: 32.67 MiB)
	# 53.091 ms (7 allocations: 48.03 MiB)

	# conclusion: dtft_loop_n is the fastest
	# (also amenable to parfor=pmap but unhelpful!?)
=#

	true
end


# 2D verify consistency
function dtft_test2c( ;
		M::Int=2^9,
		N::Dims=(2^6,2^4),
		n_shift::AbstractVector{<:Real} = zeros(Int, length(N)))

	seed!(0)
	x = randn(ComplexF64, N)
	w = (rand(M,2) .- 0.5) * 2 * pi

	o1 = dtft_loop_n(w, x; n_shift=n_shift)
	sd = dtft_init(w, N; n_shift=n_shift)
	o2 = sd.dtft(x)
	@test isequal(o2, o1)
	o3 = sd.A * x[:]
	@test isequal(o3, o1)

	b1 = dtft_adj(w, o1, N; n_shift=n_shift)
	b2 = sd.adjoint(o1)
	@test isequal(b2, b1)

	b3 = sd.A' * o1
	b3 = reshape(b3, N)
	@test isequal(b3, b1)
	true
end


"""
`dtft_test1_adj()`
test adjoint
"""
function dtft_test1_adj( ; N::Int=20, M::Int=30, n_shift::Int=5)
	seed!(0)
	w = randn(M)
	A = LinearMap{ComplexF64}(x -> dtft(w, x; n_shift=n_shift),
		y -> dtft_adj(w, y, N; n_shift=n_shift), M, N)
	isapprox(Matrix(A)', Matrix(A'))
end


"""
`dtft_test2_adj()`
test adjoint
"""
function dtft_test2_adj( ;
		M::Int=2^4,
		N::Dims=(2^3,2^2),
		n_shift::AbstractVector{<:Real} = [6,1.7])
	seed!(0)
	w = (rand(M,2) .- 0.5) * 2 * pi
	sd = dtft_init(w, N; n_shift=n_shift)
	A = sd.A
	isapprox(Matrix(A)', Matrix(A'))
end


"""
`dtft(:test)`
self test
"""
function dtft(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))

	@test dtft_test1()
	@test dtft_test1c()
	@test dtft_test1_adj()

	@test dtft_test2()
	@test dtft_test2c()
	@test dtft_test2_adj()
	true
end


# @test dtft(:test)
