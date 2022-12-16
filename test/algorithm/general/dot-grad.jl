# test/algorithm/general/dot-grad.jl

using MIRT: make_dot_gradf
using LinearAlgebra: norm, dot
using Unitful: m, s
using Test: @test, @testset, @test_throws, @inferred


@testset "dot-grad" begin
    dims = (3,4)
    Tx = Float32

    f(x) = sum(sin, x)
    g1(x) = cos.(x)
    w = Array{Tx}(undef, dims)
    function g2!(w,x)
        @. w = cos(x)
        return w
    end

    function g3(x)
        @. w = cos(x) # closure
        return w
    end
    
    x = randn(Tx, dims)
    v = randn(Tx, dims)

    @test g1(x) == g2!(w,x) == g3(x)

    dv1 = @inferred make_dot_gradf(g1)
    dv2 = @inferred make_dot_gradf(g2!, x)
    dv3 = @inferred make_dot_gradf(g3)

    @test dv1(v,x) == dv2(v,x) == dv3(v,x)
end


@testset "dot-grad-cmpl" begin
    dims = (3,4)
    Tx = ComplexF32

    wt = rand(dims...)
    f(x) = 0.5 * sum(wt .* abs2.(x))
    g1(x) = wt .* x
    w = Array{Tx}(undef, dims)
    function g2!(w, x)
        @. w = wt * x
        return w
    end

    function g3(x)
        @. w = wt * x # closure
        return w
    end
    
    x = randn(Tx, dims)
    v = randn(Tx, dims)

    @test g1(x) ≈ g2!(w,x) ≈ g3(x)

    dv1 = @inferred make_dot_gradf(g1)
    dv2 = @inferred make_dot_gradf(g2!, x)
    dv3 = @inferred make_dot_gradf(g3)

    @test dv1(v,x) ≈ dv2(v,x) ≈ dv3(v,x)
end


@testset "dot-grad-unit" begin
    dims = (3,4)
    Tx = typeof(1f0s)

    a = 1f0/s # think Hz
    b = 2f0m
    f(x) = b * sum(sin, a*x)
    Tf = typeof(oneunit(typeof(b)))
    g1(x) = a * b * cos.(a * x)

    Tw = typeof(oneunit(Tf) / oneunit(Tx))
    w = Array{Tw}(undef, dims)

    function g2!(w,x)
        @. w = a * b * cos(a * x)
        return w
    end

    function g3(x)
        @. w = a * b * cos(a * x)
        return w
    end

    x = rand(Tx, dims)
    v = rand(Tx, dims)

    @test g1(x) ≈ g2!(w,x) ≈ g3(x)

    dv1 = @inferred make_dot_gradf(g1)
    dv2 = @inferred make_dot_gradf(g2!, x, Tf)
    dv3 = @inferred make_dot_gradf(g3)

    @test dv1(v,x) ≈ dv2(v,x) ≈ dv3(v,x)
end
