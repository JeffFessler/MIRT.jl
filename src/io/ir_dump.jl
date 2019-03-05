# ir_dump.jl
# 2019-03-05 Jeff Fessler

using Printf

"""
ir_dump(x::DataType)

show all the filelds of a structure more nicely than dump() does
"""
function ir_dump(y::Any)
	x = typeof(y)
	print(x)
	fields = fieldnames(x)
	fieldtypes = x.types
	for idx in 1:length(fields)
		println()
		fd = fields[idx]
		print(stdout, " ", fd, "::")
		ft = fieldtypes[idx]
		print(ft)
		if ft <: Number
			print(stdout, " ", getfield(y, fd))
		end
		if ft == String
			print(stdout, " '", getfield(y, fd), "'")
		end
	end
	nothing
end
