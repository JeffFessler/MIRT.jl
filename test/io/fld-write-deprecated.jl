# fld-write.jl

using MIRT: fld_write
using Test: @testset, @test_throws

@testset "fld_write" begin
	@test_throws String fld_write()
end
