# ir_dump.jl

using MIRT: ir_dump

using Test: @test


x = (
    int=1, real=2., string="String",
    shorttuple=(1,2), longtuple=(1,2,3,4),
    aabool = trues(1),
)
io = isinteractive() ? stdout : IOBuffer()
@test ir_dump(io, x) isa Nothing
