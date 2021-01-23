#=
fld-write.jl
write AVS .fld to file
2019-05-12, Jeff Fessler, University of Michigan
=#

export fld_write


"""
`fld_write(file, data ; varargin)`

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

function fld_write(
	file::String,
	data::AbstractArray{<:Real} ;
	check::Bool = true,
	dir::String = "",
	endian::Symbol = :le,
	head::Array{String} = empty([""]),
#	how::Symbol = :data,
	warn::Bool = true,
	raw::Bool = false,
)

	data = fld_write_data_fix(data ; warn) # make data suitable for writing in .fld

	endian != :le && endian != :be && throw("endian '$endian' unknown")

	typedict = Dict([
		(Float32, endian === :le ? "float_le" : "xdr_float"),
		(Float64, endian === :le ? "double_le" : "xdr_double"),
		(UInt8, "byte"),
		(Int16, endian === :le ? "short_le" : "short_be"),
		(Int32, endian === :le ? "int_le" : "xdr_int"),
		])

	datatype = typedict[eltype(data)] # throws error if not there

	file = joinpath(dir, file)
#	@show file

	check && isfile(file) && throw("file '$file' exists")

	if raw # if writing data to separate ".raw" file, ensure it does not exist
		fileraw = file[1:(end-3)] * "raw"
		check && isfile(fileraw) && throw("file $fileraw exists")
	end

	# open output avs file for writing
	fid = open(file, "w")
	fraw = raw ? open(fileraw, "w") : fid

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
	host_is_le = () -> ENDIAN_BOM == 0x04030201
	fun = (host_is_le() == (endian === :le)) ?
		identity : # host/file same endian
		(host_is_le() && (endian === :be)) ?
		hton :
		(!host_is_le() && (endian === :le)) ?
		write(fraw, htol.(data)) :
		throw("not done")
	write(fraw, fun.(data))

	close(fid)

	if raw
		close(fraw)
	end

end


"""
`data = fld_write_data_fix(data)`
convert data to format suitable for writing to .fld file
"""
function fld_write_data_fix(data::AbstractArray{BigFloat} ; warn::Bool=true)
	warn && (@warn "BigFloat downgraded to Float64")
	return Float64.(data)
end

function fld_write_data_fix(data::AbstractArray{Float16} ; warn::Bool=true)
	warn && (@warn "Float16 promoted to Float32")
	return Float32.(data)
end

function fld_write_data_fix(
	data::AbstractArray{T} ;
	warn::Bool=true,
) where {T <: Union{BigInt, Int64}}
	warn && (@warn "$T downgraded to Int32")
	return Int32.(data)
end

function fld_write_data_fix(data::AbstractArray{Bool} ; warn::Bool=true)
	warn && (@warn "Bool promoted to UInt8")
	return UInt8.(data)
end

@inline fld_write_data_fix(data::AbstractArray{T} ; warn::Bool=true) where
	{T <: Union{Float32, Float64, UInt8, Int16, Int32}} = data

fld_write_data_fix(data::AbstractArray{T} ; warn::Bool=true) where {T <: Any} =
		throw(ArgumentError("unsupported type $T"))
