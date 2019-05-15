# data/z-test.jl

using Test

@test extrema(ir_load_brainweb_t1_256()) == (0.0f0, 241.0f0)
