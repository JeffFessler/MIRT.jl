# test/Afft.jl

using MIRT: Afft, embed

using LinearMapsAA: LinearMapAM, LinearMapAO
using Test: @test, @testset
using FFTW: fft, bfft

@testset "AM" begin
    samp = trues(3,2); samp[2] = false

    A = Afft(samp ; operator=false)

    @test A isa LinearMapAM

    @test A * vec(ones(3,2)) == [6, 0, 0, 0, 0]
    @test A' * [1, 0, 0, 0, 0] == vec(ones(3,2))

    @test Matrix(A)' ≈ Matrix(A')
    @test A.name == "fft"
    @test eltype(A) == ComplexF32
end


@testset "AO" begin
    samp = trues(3,2); samp[2] = false
    A = Afft(samp)
    @test A isa LinearMapAO
    @test A * ones(3,2) == [6, 0, 0, 0, 0]
    @test A' * [1, 0, 0, 0, 0] == ones(3,2)
    @test Matrix(A)' ≈ Matrix(A')
end


@testset "AO-opts" begin
    dimx = (4,3,2)
    samp = trues(dimx); samp[3] = false
    dims = (1, 3)
    T = ComplexF32
    x = randn(T, dimx)
    A = Afft(samp ; dims, fft_forward=false)
    @test A isa LinearMapAO
    @test Matrix(A)' ≈ Matrix(A')
    @test A * x ≈ bfft(x, dims)[samp] # because fft_forward=false
    x = randn(real(T), dimx)
    @test A * x ≈ bfft(x, dims)[samp]
    y = randn(T, count(samp))
    @test A' * y ≈ fft(embed(y,samp), dims)
    y = randn(real(T), count(samp))
    @test A' * y ≈ fft(embed(y,samp), dims)
    @test A.name == "fft"
    @test A.fft_forward == false
    @test A.dims == dims
    @test eltype(A) == ComplexF32
end
