# fbp/z-test.jl

using Test: @testset

list = [
"image_geom.jl"
"sino_geom.jl"
#
"cuboid_im.jl"
"disk-phantom.jl"
"ellipse_im.jl"
"ellipsoid_im.jl"
"ellipse_sino.jl"
"rect_im.jl"
"rect_sino.jl"
"rotate2d.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
