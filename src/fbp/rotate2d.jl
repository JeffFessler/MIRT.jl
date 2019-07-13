#=
rotate2d.jl
2019-03-05, Jeff Fessler, University of Michigan
=#

export rotate2d

using Test: @test


"""
`(xr,yr) = rotate2d(x, y, theta)`
2D rotation
"""
function rotate2d(x, y, theta)
	xr = cos(theta) * x + sin(theta) * y
	yr = -sin(theta) * x + cos(theta) * y
	return (xr, yr)
end


"""
`rotate2d(:test)`
self test
"""
function rotate2d(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	@test isapprox([rotate2d(1, 1, pi/2)...], [(1.,-1.)...])
	true
end


# rotate2d(:test)
