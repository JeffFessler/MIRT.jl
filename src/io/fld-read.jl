#=
fld-read.jl
Jeff Fessler
=#

export fld_header, fld_read


"""
    fld_header()

Deprecated: use `AVSfldIO.fld_header` instead.
"""
function fld_header()
    throw("fld_header is deprecated: use AVSfldIO.fld_header instead.")
end


"""
    fld_read()

Deprecated: use `FileIO.load` instead (see `AVSfldIO.fld_read`).
"""
function fld_read()
    throw("fld_read is deprecated: use FileIO.load instead.")
end
