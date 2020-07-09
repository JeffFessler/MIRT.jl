# ir_dump.jl

using MIRT: ir_dump

using Test: @test


x = (a=1, b=2)
@test ir_dump(x, io=IOBuffer()) isa Nothing
