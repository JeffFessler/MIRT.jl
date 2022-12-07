# test/fbp/rotate2d.jl

using MIRT: rotate2d
using Test: @test

@test [rotate2d(1, 1, π/2)...] ≈ [(1.,-1.)...]
