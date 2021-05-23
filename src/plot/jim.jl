#=
jim.jl
jiffy image display
2019-02-23 Jeff Fessler, University of Michigan
=#

export jim

import MIRTjim # jim

@deprecate jim(a ; kwargs...) MIRTjim.jim(a ; kwargs...)
@deprecate jim(a, b ; kwargs...) MIRTjim.jim(a, b ; kwargs...)
@deprecate jim(a, b, c ; kwargs...) MIRTjim.jim(a, b, c ; kwargs...)
@deprecate jim(a, b, c, d ; kwargs...) MIRTjim.jim(a, b, c, d ; kwargs...)
