#=
read_rdb_hdr.jl
based on
https://gitlab.com/fMRI/toppe/blob/master/+toppe/+utils/read_rdb_hdr.m
=#

export read_rdb_hdr

include("rdb-26_002.jl")

"""
`s = read_rdb_hdr(fid::IOStream)`

Read GE raw (RDB) header for MRI scans

Returns NamedTuple `s` with header values accessible by `s.key`

Matlab version is:
Copyright (c) 2012 by General Electric Company. All rights reserved.

2019-05-22 Julia version by Jeff Fessler
"""
function read_rdb_hdr(fid::IOStream)
	seek(fid, 0)
	rdbm_rev = read(fid, Float32) # raw header (RDBM) revision number

	if rdbm_rev == Float32(26.002)
		s = read_rdb_hdr_26_002(fid)
	else
		throw("unknown RDBM rev $rdbm_rev")
	end

	if (s.rdbm_rev != rdbm_rev)
		throw("rev mismatch: $(s.rdbm_rev) != $rdbm_rev")
	end

	return s
end


"""
`s = read_rdb_hdr(file::String)`
read from `file`
"""
function read_rdb_hdr(file::String)
	return open(read_rdb_hdr, file)
end


"""
`read_rdb_hdr(:test)`
self test
"""
function read_rdb_hdr(test::Symbol)
	file = "/n/ir71/d3/fessler/fmri-data-michelle-L+S/P97792.7"
	if isfile(file)
		return read_rdb_hdr(file).dab[2] == 31
	end
	true
end
