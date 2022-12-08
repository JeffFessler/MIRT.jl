# test/utility/rect.jl

using MIRT: rect

using Test: @test, @inferred

@test (@inferred rect(Int(3))) === zero(Float32)
@test (@inferred rect(Int(0))) === one(Float32)
@test (@inferred rect(0.25)) === one(Float64)
@test (@inferred rect(Float16(0.25))) === one(Float32)
