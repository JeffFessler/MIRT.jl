# downsample.jl

using MIRT: downsample_dim1
using MIRT: downsample1
using MIRT: downsample2
using MIRT: downsample3
using MIRT: downsample
import MIRT: downsample3_loop, downsample3_perm

using Test: @test, @inferred


"""
    downsample_dim1_test()
"""
function downsample_dim1_test()
	@inferred downsample_dim1(1:4, 2)
	x = reshape(2:2:48, 4, 6)
#	y = @inferred downsample_dim1(x, 2) # todo: fails
	y = downsample_dim1(x, 2)
	@test y == reshape(3:4:47, 2, 6)
	true
end


"""
downsample1_test()
"""
function downsample1_test()
	x = 2:2:48
	y = @inferred downsample1(x, 2)
	@test y == 3:4:47
	@inferred downsample1(1:5, 2, warn=false) # warn
	true
end


"""
downsample2_test()
"""
function downsample2_test()
	x = reshape(0:2:46, 4, 6)
	y = @inferred downsample2(x, 2)
	@test y == [5 21 37; 9 25 41]
	true
end


"""
downsample3_test()
"""
function downsample3_test()
	x = [6, 5, 2]
	x = reshape(2*(1:prod(x)), x...)
	for down = 1:2
		y = downsample3(x, down)
	#	y = @inferred downsample3(x, down) # todo fails
		if down == 1
			@test y == x
		elseif down == 2
			@test y == reshape([39 63; 43 67; 47 71], 3,2,1) # squeeze
		end
	end

	x = ones(4, 6, 8)
	@test downsample3_loop(x, [2, 2, 2]) == downsample3_perm(x, (2, 2, 2))

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


@test downsample_dim1_test()
@test downsample1_test()
@test downsample2_test()
@test downsample3_test()
