# nufft/z-test.jl

using Test: @test

@test dtft(:test)
@test nufft(:test)
