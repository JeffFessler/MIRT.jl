#=
reverser.jl
based on
https://stackoverflow.com/questions/27411401/julia-reverse-n-dimensional-arrays
=#

using MIRT: reverser

using Test: @test

@test reverser(1:3) == 3:-1:1
@test reverser(1:3, 1) == 3:-1:1
@test reverser((1:3)', 1) == (1:3)'
@test reverser((1:3)', 2) == (3:-1:1)'
