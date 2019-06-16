# system/z-test.jl

using Test

@test fatrix_block(:test)
@test ir_lm_test_getindex_cumsum()
@test hcat_lm(:test)
