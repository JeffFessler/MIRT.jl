# io/z-test.jl

using Test

@test ir_dump(:test)
@test fld_read(:test)
