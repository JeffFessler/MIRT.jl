# downsample.jl
# Copyright 2019-03-05, Jeff Fessler, University of Michigan

using Printf
using Test

"""
`y = downsample1(x, down; warn=true)`
downsample by factor m along first dimension by averaging

in
* x	[n1 (Nd)]
* down			integer downsampling factor

option
`warn`	warn if noninteger multiple; default true

out
* y	[n1/down (Nd)]

Copyright 2019-03-05, Jeff Fessler, University of Michigan
"""
function downsample1(x, down::Integer; warn::Bool=true)

	dim = size(x)
	dim1 = dim[1]
	x = reshape(x, dim1, :) # [n1 *Nd]
	m1 = Int(floor(dim1 / down))
	if m1 * down < dim1
		if warn
			@warn(@sprintf("truncating input size %d to %d", dim1, m1 * down))
		end
		x = x[1:(m1*down),:]
	end
	y = reshape(x, down, :)
	y = sum(y, dims=1) / down # mean
	y = reshape(y, m1, dim[2:end]...)
	return y
end


"""
downsample1_test()
"""
function downsample1_test()
	x = reshape(2:2:48, 4, 6)
	y = downsample1(x, 2)
	@test y == reshape(3:4:47, 2, 6)
	true
end


"""
`y = downsample2(x, down, warn=true)`

downsample by averaging by integer factors
in
* `x` [nx ny]
* `down` can be a scalar (same factor for both dimensions) or a 2-vector

option
* `warn`	warn if noninteger multiple; default true

out
`y`	[nx/down ny/down]
"""
function downsample2(x, down; warn::Bool=true)
	fun = (x, d) -> downsample1(x, d, warn=warn)
	y = fun(x, down[1])
	y = fun(y', down[end])'
	return y
end



"""
`downsample2_test()`
"""
function downsample2_test()
	x = reshape(1:24, 4, 6)
	y = downsample2(x, 2)
	@test y ==  [3.5 11.5 19.5; 5.5 13.5 21.5]
	true
end



"""
`y = downsample3(x, down, warn=true)`

downsample by averaging by integer factors
in
* `x` [nx ny nz]
* `down` can be a scalar (same factor for all dimensions) or a 3-vector

option
* `warn`	warn if noninteger multiple; default true

out
`y`	[nx/down ny/down nz/down]
"""
function downsample3(x, down; warn::Bool=true)

	if ndims(x) == 2
		@warn("2d case")
		if length(down) == 1; down = [down, down, 1]; end
		length(down) != 3 & throw(DimensionMismatch("$(length(down)) != 3"))&
		return downsample2(x, down[1:2])
	end

	ndims(x) != 3 && throw(DimensionMismatch("ndims(x)=$(ndims(x)) != 3"))

	if length(down) == 1
		down = [down,down,down]
	end

	length(down) != ndims(x) && throw(DimensionMismatch("length(down) != ndims(x)"))

	# downsample along each dimension
	y = downsample1(x, down[1])
	y = downsample1(permutedims(y, [2, 1, 3]), down[2])
	y = downsample1(permutedims(y, [3, 2, 1]), down[3]) # [3 1 2] order
	y = permutedims(y, [2, 3, 1])
	return y
end



function downsample3_test()
	x = [6, 5, 2]
	x = reshape(2*(1:prod(x)), x...)
	for down = 1:2
		y = downsample3(x, down)
		if down == 1
			@test y == x
		elseif down == 2
			@test y == reshape([39 63; 43 67; 47 71], 3,2,1) # squeeze
		end
	end
	true
end

#=
if has_aspire # todo
	filex = [test_dir filesep 'testx.fld']
	filey = [test_dir filesep 'testy.fld']
	fld_write(filex, x)
	delete(filey)
	com = ['op sample3 mean ' filey ' ' filex ' %d %d %d']
	com = sprintf(com, down, down, down)
	os_run(com)
	z = fld_read(filey)
	!isequal(y, z) && throw("aspire/matlab mismatch")
end
=#


"""
`downsample(:test)` self test
"""
function downsample(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	downsample1_test()
	downsample2_test()
	downsample3_test()
	true
end
