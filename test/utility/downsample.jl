# test/utility/downsample.jl

using MIRT: downsample_dim1
using MIRT: downsample1
using MIRT: downsample2
using MIRT: downsample3
import MIRT: downsample3_loop, downsample3_perm

using Test: @test, @testset, @inferred


@testset "downsample_dim1" begin
    @inferred downsample_dim1(1:4, 2)
    x = reshape(2:2:48, 4, 6)
    y = @inferred downsample_dim1(x, 2)
    y = downsample_dim1(x, 2)
    @test y == reshape(3:4:47, 2, 6)

    # stress test type inference by Int in but Float out
    x = reshape(1:6, 2, 3)
    y = @inferred downsample_dim1(x, 2)
    @test y == [1.5 3.5 5.5]
end


@testset "downsample1" begin
    x = 2:2:48
    y = @inferred downsample1(x, 2)
    @test y == 3:4:47
    @inferred downsample1(1:5, 2, warn=false) # warn

    x = 1:6
    y = @inferred downsample1(x, 2)
    @test y == 1.5:2:5.5
end


@testset "downsample2" begin
    x = reshape(0:2:46, 4, 6)
    y = @inferred downsample2(x, 2)
    @test y == [5 21 37; 9 25 41]

    x = reshape(1:12, 4, 3)
    y = @inferred downsample2(x, 2)
    @test y == reshape([3.5; 5.5], 2, 1)
end


@testset "downsample3" begin
    x = [6, 5, 2]
    x = reshape(2*(1:prod(x)), x...)
    for down = 1:2
         y = @inferred downsample3(x, down)
        if down == 1
            @test y == x
        elseif down == 2
            @test y == reshape([39 63; 43 67; 47 71], 3,2,1) # squeeze
        end
    end

    x = reshape(1:(3*4*5), 4, 3, 5)
    yl = @inferred downsample3_loop(x, (2,2,2))
    yp = @inferred downsample3_perm(x, (2,2,2))
    @test yl == yp
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
