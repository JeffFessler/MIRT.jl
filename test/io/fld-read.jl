# fld-read.jl

using MIRT: fld_header, fld_read
import MIRT: arg_get, string_to_array, fld_read_single
using Test: @test, @testset, @test_throws


dir = mktempdir()


@testset "fld eof" begin
	file = joinpath(dir, "fld-read-eof-test.fld")
	open(file, "w") do fid
		@show fid
	    println(fid, "# AVS field file eof test")
	end
	@test_throws String fld_read(file)

	open(file, "r") do fid
		@test_throws String fld_read_single(
			file, fid,
			(99,), Float32, fieldtype,
			false, # is_external_file,
			"", # extfile,
			Float32, "ieee-le", 0, 0,
		)
	end
end


@testset "fld misc" begin
	@test_throws String arg_get(["1", "2"], "3")
	@test string_to_array("0") isa AbstractArray{<:AbstractString}
end





# todo: test read slices when that feature is added
# fld_write ...
