# test/algorithm/general/ls-mm.jl

using MIRT: line_search_mm
using MIRT: make_dot_gradf, make_dot_curvf
using MIRT: LineSearchMMWork, LineSearchMMState
using LinearAlgebra: norm, dot
using Unitful: m
using Test: @test, @testset, @test_throws, @inferred
# todo: compare to LineSearches: ls HagerZhang MoreThuente BackTracking StrongWolfe Static


@testset "lsmm-work" begin
    ddims = [(3,2), (5,4)]
    uu = [randn(Float16, d) for d in ddims]
    vv = [randn(ComplexF16, d) for d in ddims]

    work = @inferred LineSearchMMWork(uu, vv, 0f0)
    @test size.(work.zz) == ddims
end


@testset "lsmm-state" begin
    ddims = [(3,2), (5,4)]
    uu = [randn(Float16, d) for d in ddims]
    vv = [randn(ComplexF16, d) for d in ddims]

    gradf = [identity, identity]
    curvf = [z -> ones(size(z)), z -> zeros(size(z))]
    @inferred make_dot_gradf(gradf[1])
    @inferred make_dot_curvf(curvf[1])
    state = @inferred LineSearchMMState(
        uu, vv, make_dot_gradf.(gradf), make_dot_curvf.(curvf),
    )
end


@testset "lsmm-iter" begin
    ddims = [(3,2), (5,4)]
    Tf = Float32
    Tc = complex(Tf)
    uu = [randn(Tf, d) for d in ddims]
    vv = [randn(Tc, d) for d in ddims]

    yy = Any[randn(Tc, ddims[1]), Tf(7)]
    aa = Tc[1+2im, 9]
    ff = [
        u -> 0.5f0 * sum(abs2, aa[1] * u .- yy[1]),
        u -> 0.5f0 * sum(abs2, aa[2] * u .- yy[2]),
    ]
    cost(α) = sum(j -> ff[j](uu[j] + α * vv[j]), eachindex(uu))
    gradf = [u -> aa[1]' * (aa[1] * u .- yy[1]),
             u -> aa[2]' * (aa[2] * u .- yy[2])]
    curvf = [z -> fill(abs2(aa[1]), size(z)),
             z -> fill(abs2(aa[2]), size(z))]
    anum = sum(j -> dot(aa[j] * vv[j], yy[j] .- aa[j] * uu[j]), 1:2)
    aden = sum(j -> norm(aa[j] * vv[j])^2, 1:2)
    ahat = real(anum) / aden

    ninner = 3 # should converge in 1 iteration since quadratic!
    out = Vector{Any}(undef, ninner+1)
    fun = (state, iter) -> (state.α)
    amm = @inferred line_search_mm(uu, vv, gradf, curvf; ninner, out, fun)
    @test amm ≈ ahat
    costs = cost.(out)
    @test cost(ahat) ≈ costs[2] ≈ costs[3] ≈ costs[4]
#   @test all(≤(0), diff(costs)) # todo
end


#=
todo: special case of a single function
=#
