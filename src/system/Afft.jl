#=
Afft.jl
2019-07-07, Jeff Fessler, University of Michigan
=#

export Afft

# using MIRT: embed
using LinearMapsAA: LinearMapAA
using FFTW: fft, ifft
using Test: @test


"""
    A = Afft(samp::AbstractArray{Bool} ; T::DataType = ComplexF32)
Make a LinearMap object for (under-sampled) FFT, of type `T`,
using given sampling pattern `samp`.
Especially for compressed sensing MRI with Cartesian sampling.
"""
function Afft(samp::AbstractArray{Bool} ; T::DataType = ComplexF32)
	dim = size(samp)
	return LinearMapAA(
		x -> fft(reshape(x,dim))[samp],
		y -> prod(dim) * vec(ifft(embed(y,samp))),
		(sum(samp), prod(dim)), (name="fft",) ; T=T,
	)
end


"""
    A = Afft(:test)
self test
"""
function Afft(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	samp = trues(3,2); samp[2] = false
	A = Afft(samp)
	@test isapprox(Matrix(A)', Matrix(A'))
	@test A * vec(ones(3,2)) == [6, 0, 0, 0, 0]
	@test A' * [1, 0, 0, 0, 0] == vec(ones(3,2))
	@test A.name == "fft"
	@test eltype(A) == ComplexF32
	true
end
