#=
fld-read.jl
Jeff Fessler
=#

export fld_header, fld_read

using AVSfldIO # temporary
import FileIO # temporary


@deprecate fld_header(args... ; kwargs...) AVSfldIO.fld_header(args... ; kwargs...)
@deprecate fld_read(args... ; kwargs...) FileIO.load(args... ; kwargs...)
