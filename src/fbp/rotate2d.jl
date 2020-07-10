#=
rotate2d.jl
2019-03-05, Jeff Fessler, University of Michigan
=#

export rotate2d


"""
`(xr,yr) = rotate2d(x, y, theta)`
2D rotation
"""
function rotate2d(x, y, theta)
	xr = cos(theta) * x + sin(theta) * y
	yr = -sin(theta) * x + cos(theta) * y
	return (xr, yr)
end
