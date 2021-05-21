#=
fld-write.jl
Jeff Fessler, University of Michigan
=#

export fld_write


"""
    fld_write()

Deprecated: use `FileIO.save` instead (see `AVSfldIO.fld_write`).
"""
function fld_write()
    throw("fld_write is deprecated: use `FileIO.save` instead")
end
