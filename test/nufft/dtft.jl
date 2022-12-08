# test/nufft/dtft.jl

import MIRT: dtft_init, dtft, dtft_adj
import MIRT: dtft_loop_n
import MIRT: dtft_loop_m
import MIRT: dtft_dist_m
import MIRT: dtft_pmap_m
import MIRT: dtft_matvec

using FFTW: fft
#using LinearAlgebra: norm
using Test: @test, @testset
using LinearMapsAA: LinearMapAA, LinearMapAO


# 1D small test data to verify correctness
function dtft_test1( ; N::Int=2^10)
    x = randn(ComplexF64, N)
    w = (2*pi) * (0:N-1) / N

    o0 = fft(x)

    o1 = dtft_loop_n(w, x)
    o2 = dtft_loop_m(w, x)
    o3 = dtft_dist_m(w, x)
    o4 = dtft_pmap_m(w, x)
    o5 = dtft_matvec(w, x)

    @test isapprox(o1, o0)
    @test isapprox(o2, o0)
    @test isapprox(o3, o0)
    @test isapprox(o4, o0)
    @test isapprox(o5, o0)

    true
end


# 2D small test
function dtft_test2( ; N::Dims = (2^3,2^2))
    x = randn(ComplexF64, N)
    w1 = (2*pi) * (0:N[1]-1) / N[1]
    w2 = (2*pi) * (0:N[2]-1) / N[2]
    w1 = repeat(w1, 1, N[2])
    w2 = repeat(w2', N[1], 1)
    w = [vec(w1) vec(w2)]
    d = dtft_init(w, N)

    o0 = fft(x)
#   o1 = dtft_loop_n(w, x)
    o1 = d.dtft(x)
    o1 = reshape(o1, N)
    @test isapprox(o1, o0)
    true
end


# 1D verify consistency
function dtft_test1c( ; N::Int=2^10, M::Int=2^11, n_shift::Real=7)
    x = randn(ComplexF64, N)
    w = randn(M)

    o1 = dtft_loop_n(w, x ; n_shift=n_shift)
    o2 = dtft_loop_m(w, x ; n_shift=n_shift)
    o3 = dtft_dist_m(w, x ; n_shift=n_shift)
    o4 = dtft_pmap_m(w, x ; n_shift=n_shift)
    o5 = dtft_matvec(w, x ; n_shift=n_shift)

    @test isapprox(o2, o1)
    @test isapprox(o3, o1)
    @test isapprox(o4, o1)
    @test isapprox(o5, o1)

    d = dtft_init(w, N ; n_shift=n_shift)
    o6 = d.dtft(x)
    @test isapprox(o6, o1)

    o7 = d.A * x
    @test isapprox(o7, o1)

    b1 = dtft_adj(w, o1, N ; n_shift=n_shift)
    b2 = d.adjoint(o1)
    @test isapprox(b2, b1)
    b3 = d.A' * o1
    @test isapprox(b3, b1)

    A = d.A
    @test A.name == "dtft1"
    @test A.N == (N,)

    true
end


# 2D verify consistency
function dtft_test2c( ;
    M::Int=2^9,
    N::Dims=(2^6,2^4),
    n_shift::AbstractVector{<:Real} = zeros(Int, length(N)),
)

    x = randn(ComplexF64, N)
    w = (rand(M,2) .- 0.5) * 2 * pi

    o1 = dtft_loop_n(w, x ; n_shift=n_shift)
    sd = dtft_init(w, N ; n_shift=n_shift)
    o2 = sd.dtft(x)
    @test isequal(o2, o1)
    o3 = sd.A * x
    @test isequal(o3, o1)

    b1 = dtft_adj(w, o1, N ; n_shift=n_shift)
    b2 = sd.adjoint(o1)
    @test isequal(b2, b1)

    b3 = sd.A' * o1
    @test isequal(b3, b1)

    A = sd.A
    @test A.name == "dtft2"
    @test A.N == N
    @test A isa LinearMapAO
    true
end


"""
    dtft_test1_adj()
test adjoint
"""
function dtft_test1_adj( ; N::Int=20, M::Int=30, n_shift::Int=5)
    w = randn(M)
    _, _, A = dtft_init(w, N ; n_shift = n_shift)
    @test Matrix(A)' ≈ Matrix(A')
    true
end


"""
    dtft_test2_adj()
test adjoint
"""
function dtft_test2_adj( ;
    M::Int=2^4,
    N::Dims=(2^3,2^2),
    n_shift::AbstractVector{<:Real} = [6,1.7],
)
    w = (rand(M,2) .- 0.5) * 2 * pi
    sd = dtft_init(w, N ; n_shift=n_shift)
    A = sd.A
    @test A isa LinearMapAO
    @test Matrix(A)' ≈ Matrix(A')
    true
end


@testset "1d" begin
    @test dtft_test1()
    @test dtft_test1c()
    @test dtft_test1_adj()
end

@testset "2d" begin
    @test dtft_test2()
    @test dtft_test2c()
    @test dtft_test2_adj()
end
