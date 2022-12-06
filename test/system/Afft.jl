# test/Afft.jl

using MIRT: Afft, embed

using LinearMapsAA: LinearMapAM, LinearMapAO
using LinearAlgebra: mul!
using Test: @test, @testset
using FFTW: fft, bfft


@testset "Afft-AM" begin # no samp
    xdim = (3,2)
    T = ComplexF32
    A = Afft(xdim ; T, operator=false)

    @test A isa LinearMapAM

    x = ones(xdim)
    y = zeros(xdim); y[1] = prod(xdim)
    @test A * vec(x) == vec(y)
    @test A' * [1, 0, 0, 0, 0, 0] == vec(ones(3,2))
    @test A.name == "fft"
    @test eltype(A) == T

    x = randn(T, xdim)
    y = fft(x)
    b = bfft(y)
    @test A * vec(x) == vec(y)
    @test A' * vec(y) == vec(b)
    @test mul!(similar(vec(y)), A, vec(x)) == vec(y)
    @test mul!(similar(vec(b)), A', vec(y)) == vec(b)

    x = randn(ComplexF64, xdim)
    y = fft(x)
    b = bfft(y)
    @test T.(A * vec(x)) ≈ T.(vec(y)) # only 32-bit accuracy!
    @test T.(A' * vec(y)) ≈ T.(vec(b))
    @test T.(mul!(similar(vec(y)), A, vec(x))) ≈ T.(vec(y))
    @test T.(mul!(similar(vec(b)), A', vec(y))) ≈ T.(vec(b))
end


@testset "Afft-AX" begin # no samp
    xdim = (3,2)
    T = ComplexF32
    fdim = 1
    sqrtN = sqrt(prod(xdim[fdim]))
    x1 = randn(T, xdim)
    T2 = ComplexF64 # for testing promotion
    x2 = randn(T2, xdim)
    for unitary in (false, true), operator in (false, true)
        A = Afft(xdim, fdim ; T, unitary, operator)
        factor = unitary ? 1/sqrtN : 1
        y1 = fft(x1, fdim) * factor
        y2 = fft(x2, fdim) * factor
        b1 = bfft(y1, fdim) * factor
        b2 = bfft(y2, fdim) * factor
        @test Matrix(A)' ≈ Matrix(A')

        if operator
            @test A * x1 ≈ y1
            @test T.(A * x2) ≈ T.(y2) # 32-bit only
        else
            @test A * vec(x1) ≈ vec(y1)
            @test T.(A * vec(x2)) ≈ T.(vec(y2)) # 32-bit only
        end
    end
end


@testset "AM" begin
    xdim = (3,2)
    samp = trues(xdim); samp[2] = false

    A = Afft(samp ; operator=false)

    @test A isa LinearMapAM

    @test A * vec(ones(xdim)) == [6, 0, 0, 0, 0]
    @test A' * [1, 0, 0, 0, 0] == vec(ones(xdim))

    @test Matrix(A)' ≈ Matrix(A')
    @test A.name == "fft"
    @test eltype(A) == ComplexF32
end


@testset "AO" begin
    xdim = (3,2)
    samp = trues(xdim); samp[2] = false
    T = ComplexF32
    A = Afft(samp; T)
    @test A isa LinearMapAO
    @test A * ones(xdim) == [6, 0, 0, 0, 0]
    @test A' * [1, 0, 0, 0, 0] == ones(3,2)
    @test Matrix(A)' ≈ Matrix(A')

    x = randn(ComplexF64, xdim)
    y = fft(x)[samp]
    b = bfft(embed(y, samp))
#   @test A * x ≈ y # nope, 32-bit only
    @test T.(A * x) ≈ T.(y)
    @test T.(A' * y) ≈ T.(b)
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
