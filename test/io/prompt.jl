# prompt.jl

import MIRT # prompt
using Test: @test

@test MIRT.prompt(:state) isa Symbol
