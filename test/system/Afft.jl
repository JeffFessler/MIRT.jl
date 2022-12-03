# test/Afft.jl

using MIRT: Afft

using LinearMapsAA: LinearMapAM, LinearMapAO
using Test: @test, @testset


samp = trues(3,2); samp[2] = false

@testset "AM" begin
	A = Afft(samp ; operator=false)
	@test A isa LinearMapAM
	@test Matrix(A)' ≈ Matrix(A')
	@test A * vec(ones(3,2)) == [6, 0, 0, 0, 0]
	@test A' * [1, 0, 0, 0, 0] == vec(ones(3,2))
	@test A.name == "fft"
	@test eltype(A) == ComplexF32
end

@testset "AO" begin
	A = Afft(samp)
	@test A isa LinearMapAO
	@test A * ones(3,2) == [6, 0, 0, 0, 0]
	@test A' * [1, 0, 0, 0, 0] == ones(3,2)
	@test Matrix(A)' ≈ Matrix(A')
end
