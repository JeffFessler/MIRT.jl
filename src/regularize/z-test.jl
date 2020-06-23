# regularize/z-test.jl

using Test: @test

# Aodwt.jl
@test Aodwt(:test)

# diffl.jl
@test diffl(:test)
@test diffl_map(:test)

# diffs.jl
@test diff_map(:test)
