#=
fld-write.jl
Jeff Fessler, University of Michigan
=#

export fld_write

import FileIO # temporary

@deprecate fld_write(args... ; kwargs...) FileIO.save(args... ; kwargs...)
