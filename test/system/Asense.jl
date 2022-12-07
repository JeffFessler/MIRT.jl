# test/system/Asense.jl

using MIRT: Asense, embed
using LinearMapsAA: LinearMapAO
using FFTW: fft, bfft, fftshift, ifftshift
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
    smaps = randn(T, dims..., ncoil)
    A = Asense(samp, smaps)
    x = randn(T, dims)
    y = randn(T, count(samp), ncoil)
    @test isapprox(dot(y, A * x), dot(A' * y, x))
end

@testset "Asense-opt" begin
end
    sdim = (3,4,2)
    samp = trues(sdim); samp[2] = false
    T = ComplexF32
    smaps = [randn(T, sdim), randn(T, sdim)]
    dims = [1, 3]
    sqrtN = sqrt(prod(sdim[dims]))
    for unitary in (false, true)
        A = Asense(samp, smaps; dims, fft_forward=false, unitary)
        factor = unitary ? 1/sqrtN : 1
        @test Matrix(A)' ≈ Matrix(A')

        # note bfft here because fft_forward = false, for code coverage
        forw(x, smap) = fftshift(bfft(ifftshift(x .* smap), dims))[samp] * T(factor)
        x = randn(T, sdim)
        y = [forw(x, smaps[1]) forw(x, smaps[2])]
        @test A * x ≈ y

        back(z, smap) = fftshift(fft(ifftshift(z), dims)) .* conj(smap) * T(factor)
        b = back(embed(y[:,1], samp), smaps[1]) +
            back(embed(y[:,2], samp), smaps[2])
        @test A' * y ≈ b
    end
