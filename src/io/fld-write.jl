# fld-write.jl
# write AVD .fld to file
# 2019-05-12, Jeff Fessler, University of Michigan

"""
`dir = ir_test_dir!(path)`
set default test directory
"""
function ir_test_dir!(path::String)
	ENV["TEST_DIR"] = path
end


"""
`dir = ir_test_dir`
returns default test directory; default is `/tmp`
"""
function ir_test_dir()
#	get(ENV, "TEST_DIR", "/tmp")
	try
		return ENV["TEST_DIR"]
	catch
		return "/tmp"
	end
end


"""
`fld_write(file, data, varargin)`

write data into AVS format `.fld` file

See ASPIRE user's guide (or AVS online documentation) for file format.

in
* `file` name of file typically ending in `.fld`
* `data` real data array

option
* `check`	Bool	report error if file exists; default `true`
* `dir`		String	directory name to prepend file name; default `""`
* `head`	Array{String}	comment information for file header
* `how`		Symbol	what file data type; default `data` infer from data
					todo: e.g., 'short_be'.  default is 'xdr_float'
* `raw`		Bool	put raw data in `name.raw`, header in `name.fld`
					where `file` = `name.fld`; default `false`

At this point, only little endian is supported
* `endian`	String	'ieee-le' or 'ieee-be' for little/big endian
				default inferred from 'how'
				with bias to 'ieee-le' - todo: not yet implemented
"""

function fld_write(file::String, data::AbstractArray{<:Real};
	check::Number = true,
	dir::String = "",
	how::Symbol = :data,
	head::Array{String} = empty([""]),
	raw::Bool = false,
#	endian:: = ?, # infer below from arg.type
	)


#if isempty(arg.endian)
#	switch arg.type
#case {'short_be', 'int_be', 'float_be', 'double_be', ...
#	'xdr_int', 'xdr_double' 'xdr_float'}
#		arg.endian = 'ieee-be'
#	otherwise
#		arg.endian = 'ieee-le'
#	end
#end

	data = fld_write_data_fix(data) # make data suitable for writing in .fld

	typedict = Dict([
		(Float32, "float_le"),
		(Float64, "double_le"),
		(Int8, "byte"),
		(Int16, "short_le"),
		(Int32, "int_le"),
		])

	dtype = eltype(data)
	datatype = typedict[dtype] # throws error if not there

	file = joinpath(dir, file)
#	printm('file = "%s"', file)

	check && isfile(file) && throw("file $file exists")

	# write data to separate ".raw" file?
	if raw
		fileraw = file[1:(end-3)] * "raw"
		check && isfile(fileraw) && throw("file $fileraw exists")
	end

	# open avs file for writing
	fid = open(file, "w")
	fraw = fid # default
	if raw
        fraw = open(fileraw, "w")
	end

	# write header
	ndim = ndims(data)

	println(fid, "# created by $(basename(@__FILE__))")
	for ii=1:length(head)
		println(fid, "# $(head[ii])")
	end

	println(fid, "ndim=$ndim")
#	println(fid, "nspace=$ndim")
	for ii=1:ndim
		println(fid, "dim$ii=$(size(data,ii))")
	end
	println(fid, "data=$datatype")
	println(fid, "veclen=1")
	println(fid, "field=uniform")

	if raw
		println(fid, "variable 1 file=$fileraw filetype=binary")
	else
		write(fid, "\f\f") # two form feeds: char(12)
	end

	# finally, write the binary data
	write(fraw, data)

	close(fid)

	if raw
		close(fraw)
	end

end


"""
`data = fld_write_data_fix(data)`
convert data to format suitable for writing to .fld file
"""
function fld_write_data_fix(data::AbstractArray{<:Real})

	dtype = eltype(data)

	dtype <: Unsigned && throw("Only Real and Signed types supported")

	if dtype == BigFloat
		@warn "BigFloat downgraded to Float64"
		return Float64.(data)
	end

	if dtype == Float16
		@warn "Float16 upgraded to Float32"
		return Float32.(data)
	end

	if dtype == BigInt
		@warn "BigInt downgraded to Int32"
		return Int32.(data)
	end

	if dtype == Int128
		@warn "Int128 downgraded to Int64"
		return Int64.(data)
	end

	if dtype == Int64
		@warn "Int64 downgraded to Int32"
		return Int32.(data)
	end

	if eltype(data) == Bool
		@warn "Bool upgraded to Int8"
		return Int8.(data)
	end

	return data
end



"""
`fld_write_test1(file, data, ...)`
"""
function fld_write_test1(file, data; raw::Bool=false, chat::Bool=false, kwarg...)
	fld_write(file, data; raw=raw, kwarg...)
	tmp = fld_read(file, chat=chat)
	tmp != data && throw("test failed for file = $file")

	rm(file)
	if raw
		tmp = file[1:(end-4)] * ".raw"
		rm(tmp)
	end
	true
end


"""
`fld_write(:test)`
"""
function fld_write(test::Symbol; chat::Bool=false)
	test != :test && throw("bug")

	file = joinpath(ir_test_dir(), "fld-write-test.fld")

	if true
		chat && @info "test1"
		data = Float32.([5:8; 1:4])
		fld_write_test1(file, data, check=true)
	end

	if true
		chat && @info "test2"
		data = Int32.([5:8; 1:4])
		fld_write_test1(file, data, check=true)
	end

	if true # test raw/header
		chat && @info "test3"
		fld_write_test1(file, data, check=true, raw=true, chat=chat)
	end

#	formats = {'float_be', 'float_le', 'float', 'xdr_float'}
#	for ii=1:length(formats)
	#	format = formats{ii}
	#	pr format
	#	fld_write_test1(file, data, 'check', 1, 'type', format)
	#	%		 'endian', 'ieee-le', ... % nah, defer from type
#	end
	#	format = 'short_le'
	#	pr format
	#	fld_write_test1(file, int16(data), 'check', 1, 'type', format)

	true
end
