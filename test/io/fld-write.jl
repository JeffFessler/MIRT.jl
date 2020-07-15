# fld-write.jl

using MIRT: fld_write, fld_read

#using MIRT: ir_test_dir, ir_test_dir!
using Test: @test, @inferred


#=
"""
`dir = ir_test_dir!(path)`

set default test directory
"""
function ir_test_dir!(path::String)
	ENV["MIRT_TEST_DIR"] = path
end


"""
`dir = ir_test_dir`

return default test directory `ENV["MIRT_TEST_DIR"]`; default is `/tmp`
"""
function ir_test_dir()
#	get(ENV, "MIRT_TEST_DIR", "/tmp")
	try
		return ENV["MIRT_TEST_DIR"]
	catch
		return "/tmp"
	end
end
=#



"""
`fld_write_test1(file, data, ...)`
"""
function fld_write_test1(file, data; raw::Bool=false, chat::Bool=false, kwarg...)
	@inferred fld_write(file, data; raw=raw, kwarg...)
	tmp = fld_read(file, chat=chat)
	tmp != data && throw("test failed for file = $file")

	rm(file)
	if raw
		tmp = file[1:(end-4)] * ".raw"
		rm(tmp)
	end
	true
end


# test writing and reading AVS `.fld` files of all key types

dir = mktempdir()
file = joinpath(dir, "fld-write-test.fld")
chat = isinteractive()

for dtype in (UInt8, Int16, Int32, Float32, Float64)
	for endian in (:le, :be)
		for raw in (false, true)
			chat && @info "dtype=$dtype endian=$endian raw=$raw"
			data = convert.(dtype, [5:8; 1:4])
			@test fld_write_test1(file, data, endian=endian, raw=raw,
				check=true, chat=false)
		end
	end
end
