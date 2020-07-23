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


@testset "fld ascii" begin
	tdir = mktempdir()
	file = joinpath(tdir, "fld-read-test.fld")
	extfilename = "fld-read-test.txt"
	chat = isinteractive()

	data = Float32.(1:5)
	head = [
		"# AVS field file"
		"ndim=1"
		"dim1=5"
		"nspace=1"
		"veclen=1"
		"data=float"
		"field=uniform"
		"variable 1 file=$extfilename filetype=ascii"
	]
	open(file, "w") do fid
		for line in head
			println(fid, line)
		end
	end

	@test fld_header(file) isa Array{<:AbstractString}

	# ensure that it throws if no extfile
	@test_throws String tmp = fld_read(file)

	# write ascii and test
	extfile = joinpath(tdir, extfilename)
	open(extfile, "w") do fid
		for d in data
			println(fid, d)
		end
	end

#	run(`op range $file`)
	tmp = fld_read(file ; chat=chat)
	@test tmp == data
	@test eltype(tmp) == Float32

	rm(file)
	rm(extfile)
end


@testset "read no ext" begin # missing extfile
	tdir = mktempdir()
	file = joinpath(tdir, "fld-read-test.fld")
	head = [
		"# AVS field file"
		"ndim=1"
		"dim1=4"
		"nspace=1"
		"veclen=1"
		"data=float"
		"field=uniform"
		"variable 1 file=bad.txt filetype=binary"
	]
	open(file, "w") do fid
		for line in head
			println(fid, line)
		end
	end

	@test_throws String fld_read(file)
	rm(file)
end


# todo: test read slices when that feature is added
# fld_write ...
