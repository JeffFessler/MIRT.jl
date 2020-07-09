#=
fld-write.jl
write AVD .fld to file
2019-05-12, Jeff Fessler, University of Michigan
=#

export fld_write


"""
`fld_write(file, data, varargin)`

write data into AVS format `.fld` file

See ASPIRE user's guide (or AVS online documentation) for file format.

in
- `file` name of file typically ending in `.fld`
- `data` real data array

option
- `check::Bool`			report error if file exists; default `true`
- `dir::String`			directory name to prepend file name; default `""`
- `endian::`Symbol`		`:le` little endian (default), `:be` big endian
- `head::Array{String}`	comment information for file header
- `raw::Bool`			put raw data in `name.raw`, header in `name.fld`
						where `file` = `name.fld`; default `false`
"""

# * `how::Symbol`			file data type; default `:data` infer from data

function fld_write(file::String, data::AbstractArray{<:Real};
	check::Bool = true,
	dir::String = "",
	endian::Symbol = :le,
	head::Array{String} = empty([""]),
#	how::Symbol = :data,
	raw::Bool = false,
	)

	data = fld_write_data_fix(data) # make data suitable for writing in .fld

	endian != :le && endian != :be && throw("endian '$endian' unknown")

	typedict = Dict([
		(Float32, endian === :le ? "float_le" : "xdr_float"),
		(Float64, endian === :le ? "double_le" : "xdr_double"),
		(UInt8, "byte"),
		(Int16, endian === :le ? "short_le" : "short_be"),
		(Int32, endian === :le ? "int_le" : "xdr_int"),
		])

	dtype = eltype(data)
	datatype = typedict[dtype] # throws error if not there

	file = joinpath(dir, file)
#	@show file

	check && isfile(file) && throw("file '$file' exists")

	if raw # if writing data to separate ".raw" file, ensure it does not exist
		fileraw = file[1:(end-3)] * "raw"
		check && isfile(fileraw) && throw("file $fileraw exists")
	end

	# open output avs file for writing
	fid = open(file, "w")
	fraw = fid # default
	if raw
        fraw = open(fileraw, "w")
	end

	# write header
	ndim = ndims(data)

	println(fid, "# AVS field file ($(basename(@__FILE__)))")
	for line in head
		println(fid, "# $line")
	end

	println(fid, "ndim=$ndim")
	for ii=1:ndim
		println(fid, "dim$ii=$(size(data,ii))")
	end
	println(fid, "nspace=$ndim")
	println(fid, "veclen=1")
	println(fid, "data=$datatype")
	println(fid, "field=uniform")

	if raw
		println(fid, "variable 1 file=$fileraw filetype=binary")
	else
		write(fid, "\f\f") # two form feeds: char(12)
	end

	# finally, write the binary data
	host_is_le = ENDIAN_BOM == 0x04030201
	if host_is_le == (endian === :le) # host/file same endian
		write(fraw, data) # finally, write the binary data
	elseif host_is_le && (endian === :be)
		write(fraw, hton.(data))
	elseif !host_is_le && (endian === :le)
		write(fraw, htol.(data))
	end

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
		@warn "Bool upgraded to UInt8"
		return UInt8.(data)
	end

	!in(dtype, (Float32, Float64, UInt8, Int16, Int32)) &&
		throw("unsupported type $dtype")

	return data
end
