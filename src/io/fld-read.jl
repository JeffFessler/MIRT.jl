#=
fld-read.jl
Jeff Fessler and David Hong
=#

export fld_header, fld_read


"""
`head = fld_header(file::String, ...)`

`head, is_external_file, fid = fld_header(file, keepopen=true)`

read header data from AVS format `.fld` file

in
- `file::String`	file name, usually ending in `.fld`

option
- `dir::String`		prepend file name with this directory; default ""
- `chat::Bool`		verbose?
- `keepopen::Bool`	true to leave file open for more reading; default false to close

out
- `head::String`	array of header information

"""
function fld_header(file::AbstractString ;
		dir::String = "",
		chat::Bool = false,
		keepopen::Bool = false,
	)

	file = joinpath(dir, file)
	fid = open(file)

	formfeed = '\f' # form feed
	newline = '\n'
	is_external_file = false

	# read header until we find the end of file or the 1st form feed
	header = ""
	while true
		# end of file means external file (or error)
		if eof(fid)
			if occursin("file=", header)
				is_external_file = true
				break # end of header file!
			else
				chat && @info(header)
				throw("end of file before form feeds?")
			end
		end

		inchar = read(fid, Char) # read one character

		# form feed means embedded file
		if inchar == formfeed
			eof(fid) && throw("end of file before 2nd form feed?")
			inchar = read(fid, Char)
			(inchar != formfeed) && throw("not two form feeds?")
			chat && @info("embedded data file")
			break
		end

		# otherwise append this character to header string
		header *= inchar
	end

	header = string_to_array(header) # convert to array
	chat && @info(header)

	if keepopen
		return header, is_external_file, fid
	end

	# else
	close(fid)
	return header
end


# todo:
# + [ ] 'raw'	0|1	1: return raw data class (default), 0: return doubles
# + [ ] 'slice' int	specify which slice to read from 3D file. (0 ... nz-1)
# + [ ] 'chat'	0|1	enable verbosity
# + [ ] 'dim_only' 0|1	returns dims. data and coord are equal to [].
# + [ ] 'coord' 0|1	returns coordinates too (default: 0)
# + [ ] 'coord_format'	default: 'n'; see fopen() for machine format options
#	(needed for some UM .fld files with 'vaxd' coordinates)
# + [ ] multi files
# + [ ] short datatype

"""
`fld_read(file::String)`

read data from AVS format `.fld` file

in
- `file`	file name, usually ending in `.fld`

option
- `dir`		String	prepend file name with this directory; default ""
- `chat`	Bool	verbose?

out
- `data`	Array (1D - 5D) in the data type of the file itself

"""
function fld_read(file::AbstractString ;
		dir::String = "",
		chat::Bool = false,
	)

	file = joinpath(dir, file)

	header, is_external_file, fid = fld_header(file, keepopen=true, chat=chat)

	# parse header to determine data dimensions and type
	ndim = arg_get(header, "ndim")
	dims = [arg_get(header, "dim$ii") for ii in 1:ndim]
	fieldtype = arg_get(header, "field", false)
	datatype = arg_get(header, "data", false)
	(arg_get(header, "veclen") != 1) && throw("only veclen=1 done")
	chat && @info("ndim=$ndim")
	chat && @info("dims=$dims")

	# external file (binary data in another file)
	# todo: external ASCII files to be implemented (from fld_read.m)
	_skip = 0
	if is_external_file
		close(fid)
		extfile = arg_get(header, "file", false)

		filetype = arg_get(header, "filetype", false)
		chat && @info("Current file = '$file', External file = '$extfile', type='$filetype'")

		if occursin("skip=",prod(header))
			_skip = arg_get(prod(header),"skip")
		end

		if filetype != "multi"
			if !isfile(extfile)
				fdir = file
				slash = findlast(isequal('/'),fdir)
				isnothing(slash) && throw("cannot find external file $extfile")
				fdir = fdir[1:slash]
				extfile = fdir*extfile # add directory
				!isfile(extfile) && throw("no external ref file $extfile")
			end
		else
			throw("multi not supported yet")
		end
	else
		filetype = ""
		extfile = ""
	end

	# finally, read the binary data
	format, endian, bytes = datatype_fld_to_mat(datatype)

	# single file reading
	data = fld_read_single(file, fid, dims, datatype, fieldtype,
		is_external_file, extfile, format, endian, bytes, _skip)

	close(fid)

	return data
end


function fld_read_single(file, fid, dims, datatype, fieldtype,
	is_external_file, extfile, format, endian, bytes, _skip)

	# reopen file to same position, with appropriate endian too.
	if is_external_file
		fid = open(extfile)
	end

	skip(fid,_skip)

	rdims = dims # from handling slice option

	# read binary data and reshape appropriately
	data = Array{format}(undef,rdims...)
	try
		read!(fid,data)
	catch
		@info("rdims=$rdims")
		throw("file count != data count")
	end
	if endian == "ieee-le"
		data .= htol.(data)
	elseif endian == "ieee-be"
		data .= hton.(data)
	end

	return data
end


"""
`header = string_to_array(header_lines)`

convert long string with embedded newlines into string array
"""
function string_to_array(header_lines::String)
	newline = '\n'

	# ensure there is a newline at end, since dumb editors can forget...
	if header_lines[end] != newline
		header_lines *= newline
	end

	ends = findall(isequal(newline), header_lines)
	(length(ends) <= 0) && throw("no newlines?")

	header = split(header_lines, newline, keepempty=false)

	# strip comments (lines that begin with #)
	header = filter(line -> line[1] != '#', header)

	return header
end


"""
`arg_get(head, name, toint)`

Parse an argument from header, of the name=value form
"""
function arg_get(head::Array{<:AbstractString}, name::String, toint::Bool=true)
	for ll = 1:length(head)
		line = head[ll]
		start = findfirst(name * '=',line)
		if !isnothing(start)
			!isnothing(findnext(name*'=',line,start[end]+1)) && throw("bug: multiples?")
			line = line[(start[end]+1):end]
			arg = split(line)[1]
			toint && (arg = parse(Int,arg))
			return arg
		end
	end
	throw("could not find $name in header")
end


"""
`format, endian, bytes = datatype_fld_to_mat(datatype)`

determine data format from .fld header datatype
"""
function datatype_fld_to_mat(datatype::AbstractString)
	if datatype == "byte"
		format = UInt8
		endian = "ieee-be" # irrelevant
		bytes = 1

	elseif datatype in ["short_be", "short_sun", "xdr_short"]
		format = Int16
		endian = "ieee-be"
		bytes = 2
	elseif datatype == "short_le"
		format = Int16
		endian = "ieee-le"
		bytes = 2

	elseif datatype == "int"
		format = Int32
		endian = "" # native int - not portable
		bytes = 4
	elseif datatype == "int_le"
		format = Int32
		endian = "ieee-le"
		bytes = 4
	elseif datatype in ["int_be", "xdr_int"]
		format = Int32
		endian = "ieee-be"
		bytes = 4

	elseif datatype == "float"
		format = Float32 # typeof(1.)
		endian = "" # native float - not portable
		bytes = 4
	elseif datatype == "float_le"
		format = Float32
		endian = "ieee-le"
		bytes = 4
	elseif datatype in ["float_be", "xdr_float"]
		format = Float32
		endian = "ieee-be"
		bytes = 4

	elseif datatype == "double"
		format = Float64
		endian = "" # native double - not portable
		bytes = 8
	elseif datatype == "double_le"
		format = Float64
		endian = "ieee-le"
		bytes = 8
	elseif datatype in ["double_be", "xdr_double"]
		format = Float64
		endian = "ieee-be"
		bytes = 8

	else
		throw("format '$datatype' not yet implemented. ask jeff!")
	end

	return format, endian, bytes
end


"""
`fld_read(:test)`
self test
"""
function fld_read(test::Symbol)
	# todo: test read slices when that feature is added
	return fld_write(test)
end
