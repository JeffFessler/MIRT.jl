#=
jim.jl
jiffy image display
2019-02-23 Jeff Fessler, University of Michigan
=#

export jim

import MIRTjim # jim

@deprecate jim(args... ; kwargs...) MIRTjim.jim(args... ; kwargs...)
