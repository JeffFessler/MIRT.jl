#=
ir_dump.jl
2019-03-05 Jeff Fessler
=#

export ir_dump


"""
    ir_dump(x::Any ; io::IO = stdout)
    ir_dump(io::IO, x::Any)

Show all the fields of a structure or `NamedTuple` more nicely than dump() does
"""
function ir_dump(y::Any ; io::IO = stdout)
	x = typeof(y)
	print(io, x)
	fields = fieldnames(x)
	fieldtypes = x.types
	for (idx,fd) in enumerate(fields)
		println(io)
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


ir_dump(io::IO, x::Any) = ir_dump(x ; io=io)
