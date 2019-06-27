# io/z-test.jl

using Test

@test ir_dump(:test)
@test fld_read(:test)
@test fld_write(:test)
#@test loadpfile(:test)
#@test prompt(:test)
