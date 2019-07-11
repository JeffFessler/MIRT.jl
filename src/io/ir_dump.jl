#=
ir_dump.jl
2019-03-05 Jeff Fessler
=#

export ir_dump


"""
`ir_dump(x::Any ; io::IO = stdout)`

Show all the fields of a structure or `NamedTuple` more nicely than dump() does
"""
function ir_dump(y::Any ; io::IO = stdout)
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
`ir_dump(:test; io::IO = IOBuffer())`

self test
"""
function ir_dump(test::Symbol ; io=IOBuffer())
	test != :test && throw(ArgumentError("test $test"))
	x = (a=1, b=2)
	ir_dump(x, io=io)
	true
end
