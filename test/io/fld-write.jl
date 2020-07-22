# fld-write.jl

using MIRT: fld_write, fld_read
#import MIRT: fld_write_data_fix

#using MIRT: ir_test_dir, ir_test_dir!
using Test: @test, @testset, @test_throws, @inferred


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
function fld_write_test1(file, data ; raw::Bool=false, chat::Bool=false, kwarg...)
	@inferred fld_write(file, data ; raw=raw, kwarg...)
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

@testset "write std" begin
	for dtype in (UInt8, Int16, Int32, Float32, Float64)
		for endian in (:le, :be)
			for raw in (false, true)
				chat && @info "dtype=$dtype endian=$endian raw=$raw"
				data = convert.(dtype, [5:8; 1:4])
				@test fld_write_test1(file, data, endian=endian, raw=raw,
					head = ["dtype=$dtype endian=$endian raw=$raw",],
					check=true, chat=false)
			end
		end
	end
end

@testset "write odd" begin
	for T in (BigFloat, Float16, BigInt, Int64, Bool)
		data = T.([1 0 1])
		@test fld_write_test1(file, data ; check=true, chat=false)
	end
end

@testset "fld convert" begin
#=
	pairs = (
		(BigFloat, Float64),
		(Float16, Float32),
		(BigInt, Int32),
		(Int64, Int32),
		(Bool, UInt8),
	)
	for (in, out) in pairs
		data = in.([1 0 1])
		@test eltype(fld_write_data_fix(data)) == out
	end
=#
	@test_throws ArgumentError fld_write_data_fix(Rational.([1 0 1]))
end
