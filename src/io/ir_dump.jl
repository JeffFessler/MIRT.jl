# ir_dump.jl
# 2019-03-05 Jeff Fessler

using Printf



"""
`ir_dump([io::IO,] x::DataType)`

show all the filelds of a structure more nicely than dump() does
"""
function ir_dump(y::Any)
	ir_dump(stdout, y)
end

function ir_dump(io::IO, y::Any)
	x = typeof(y)
	print(io, x)
	fields = fieldnames(x)
	fieldtypes = x.types
	for idx in 1:length(fields)
		println(io)
		fd = fields[idx]
		print(io, " ", fd, "::")
		ft = fieldtypes[idx]
		print(io, ft)
		if ft <: Number
			print(io, " ", getfield(y, fd))
		end
		if ft == String
			print(io, " '", getfield(y, fd), "'")
		end
	end
	println(io)
	nothing
end


"""
`ir_dump([io::IO,] :test)`

self test
"""
function ir_dump(test::Symbol)
	ir_dump(IOBuffer(), test)
end

function ir_dump(io::IO, test::Symbol)
	@assert test == :test
	x = (a=1, b=2)
	ir_dump(io, x)
	true
end
