# nufft/z-test.jl

using Test

@test dtft(:test)
@test nufft(:test)
