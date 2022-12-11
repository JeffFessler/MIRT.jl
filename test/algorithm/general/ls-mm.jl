# test/algorithm/general/ls-mm.jl

using MIRT: line_search_mm
using MIRT: LineSearchMMWork, LineSearchMMState, _dot_gradf, _dot_curvf
using Test: @test, @testset, @test_throws, @inferred
using Unitful: m
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
    @inferred _dot_gradf(gradf[1])
    @inferred _dot_curvf(curvf[1])
    state = @inferred LineSearchMMState(
        uu, vv, _dot_gradf.(gradf), _dot_curvf.(curvf),
    )
end

@testset "lsmm-iter" begin
end
    ddims = [(3,2), (5,4)]
    uu = [randn(Float16, d) for d in ddims]
    vv = [randn(ComplexF16, d) for d in ddims]

    ff = [
        u -> 0.5 * sum(abs2, u .- 6),
        u -> 1.0 * sum(abs2, u .- 7),
    ]
    cost(α) = sum(j -> ff[j](uu[j] + α * vv[j]), eachindex(uu))

    gradf = [u -> u .- 6, u -> 2 * (u .- 7)]
    curvf = [z -> ones(size(z)), z -> 2*zeros(size(z))]

    ninner = 20
    out = Vector{Any}(undef, ninner+1)
    fun = (state, iter) -> (state.α)
    @inferred line_search_mm(uu, vv, gradf, curvf; ninner, out, fun)
    cost.(out)
#   @test all(≤(0), diff(cost)) # todo


#=
todo: special case of a single function
=#
