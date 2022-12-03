# test/Asense.jl

using MIRT: Asense
using LinearMapsAA: LinearMapAO
using LinearAlgebra: dot
using Test: @test, @testset

@testset "Asense" begin
    dims = (6,5)
    samp = trues(dims); samp[2] = false
    T = ComplexF32
    smaps = [randn(T, dims), randn(T, dims)]

	A = Asense(samp, smaps)
	@test A isa LinearMapAO
	@test Matrix(A)' ≈ Matrix(A')
	@test A.name == "Asense"
	@test eltype(A) == ComplexF32

    dims = (30,20)
    samp = rand(dims...) .< 0.5
    T = ComplexF32
    ncoil = 2
    smaps = [randn(T, dims), randn(T, dims)]
	A = Asense(samp, smaps)
    x = randn(T, dims)
    y = randn(T, count(samp), ncoil)
    @test isapprox(dot(y, A * x), dot(A' * y, x))
end
