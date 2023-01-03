# test/algorithm/general/dot-curv.jl

using MIRT: make_dot_curvf
using LinearAlgebra: norm, dot
using Unitful: m, s
using Test: @test, @testset, @test_throws, @inferred


@testset "dot-curv-cmpl" begin
    dims = (3,4)
    Tx = ComplexF32

    wt = rand(dims...)
    f(x) = 0.5 * sum(wt .* abs2.(x))
    c1(x) = wt # Hessian's diagonal
    Tw = real(Tx)
    w = Array{Tw}(undef, dims)
    function c2!(w,x)
        @. w = wt
        return w
    end

    function c3(x)
        @. w = wt
        return w
    end
    
    x = randn(Tx, dims)
    v = randn(Tx, dims)

    @test c1(x) ≈ c2!(w,x) ≈ c3(x)

    dc1 = @inferred make_dot_curvf(c1)
    dc2 = @inferred make_dot_curvf(c2!, x)
    dc3 = @inferred make_dot_curvf(c3)

    @test dc1(v,x) ≈ dc2(v,x) ≈ dc3(v,x)

    c4 = maximum(wt) # max curvature
    dc4 = @inferred make_dot_curvf(c4)
    @test dc4(v,x) ≈ c4 * sum(abs2, v)
end


@testset "dot-curv-unit" begin
    dims = (3,4)
    Tx = typeof(1f0s)

    a = rand(Float32, dims...) / s # think Hz
    b = 2f0m
    f(x) = b * sum(sin, a .* x)
    Tf = typeof(oneunit(typeof(b)))
    c1(x) = a.^2 * b .* sin.(a .* x)

    Tw = typeof(oneunit(Tf) / oneunit(Tx)^2)
    w = Array{Tw}(undef, dims)

    function c2!(w,x)
        @. w = a^2 * b * sin(a * x)
        return w
    end

    function c3(x)
        @. w = a^2 * b * sin(a * x)
        return w
    end

    x = rand(Tx, dims)
    v = rand(Tx, dims)

    @test c1(x) ≈ c2!(w,x) ≈ c3(x)

    dc1 = @inferred make_dot_curvf(c1)
    dc2 = @inferred make_dot_curvf(c2!, x, Tf)
    dc3 = @inferred make_dot_curvf(c3)

    @test dc1(v,x) ≈ dc2(v,x) ≈ dc3(v,x)

    c4 = b * maximum(abs2, a) # max curvature
    dc4 = @inferred make_dot_curvf(c4)
    @test dc4(v,x) ≈ c4 * sum(abs2, v)
end
