# prompt.jl

using MIRT: prompt
using Test: @test


_tmp = prompt(:state) # save current state
prompt(:draw)
@test prompt(:state) === :draw
prompt()
@test prompt(_tmp) isa Symbol # return to original state
