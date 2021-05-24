# jim.jl

import MIRT # jim
using Test: @test, @test_throws

@test MIRT.jim(:defs) isa Dict
