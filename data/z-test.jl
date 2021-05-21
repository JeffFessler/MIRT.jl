# data/z-test.jl

using MIRT: ir_load_brainweb_t1_256
using Test: @test, @testset

@testset brainweb begin
    x = ir_load_brainweb_t1_256()
    @test extrema(x) == (0.0f0, 241.0f0)
end
