# fld-read.jl

using MIRT: fld_header, fld_read
using Test: @testset, @test_throws

@testset "fld-read" begin
	@test_throws String fld_header()
	@test_throws String fld_read()
end
