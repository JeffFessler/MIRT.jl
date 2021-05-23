#=
prompt.jl
prompt user to hit a key to continue
=#

export prompt

import MIRTjim # prompt

@deprecate prompt(args... ; kwargs...) MIRTjim.prompt(args... ; kwargs...)
