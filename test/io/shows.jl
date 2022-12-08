# test/io/shows.jl

using MIRT: @shows
using Test: @test

var = ones(3,4)
@test (@shows var) isa Nothing
