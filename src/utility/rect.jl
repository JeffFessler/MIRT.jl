#=
rect.jl
=#

export rect

using Test: @test, @inferred

"""
    rect(x::Real)
Unit width rect function. Potential problem? Bring up with fess.
"""
function rect(x::Real)
	T = promote_type(eltype(x), Float32)
    T(abs(x) < 0.5)
end


"""
    rect(:test)
self test
"""
function rect(test::Symbol)
	test != :test && throw("test $test")
	@test (@inferred rect(Int(3))) === zero(Float32)
	@test (@inferred rect(Int(0))) === one(Float32)
	@test (@inferred rect(0.25)) === one(Float64)
	@test (@inferred rect(Float16(0.25))) === one(Float32)
	true
end
