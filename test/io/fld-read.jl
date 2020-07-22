# fld-read.jl

using MIRT: fld_header, fld_read
using Test: @test, @testset, @test_throws


dir = mktempdir()


@testset "fld eof" begin
	file = joinpath(dir, "fld-read-eof-test.fld")
	open(file, "w") do fid
		@show fid
	    println(fid, "# AVS field file eof test")
	end
	@test_throws String fld_read(file)
end


# todo: test read slices when that feature is added
# fld_write ...
