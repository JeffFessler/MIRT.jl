# fbp/z-test.jl

using Test: @testset, @test_throws
import MIRT

list = [
"rotate2d.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end

@testset "deprecate" begin
 @test_throws String MIRT.image_geom()
 @test_throws String MIRT.image_geom_mri()
 @test_throws String MIRT.cuboid_im()
 @test_throws String MIRT.disk_phantom_params()
 @test_throws String MIRT.ellipse_im()
 @test_throws String MIRT.ellipse_im_params()
 @test_throws String MIRT.ellipse_sino()
 @test_throws String MIRT.ellipsoid_im()
 @test_throws String MIRT.rect_im()
 @test_throws String MIRT.rect_sino()
 @test_throws String MIRT.sino_geom()
 @test_throws String MIRT.sino_plot()
end
