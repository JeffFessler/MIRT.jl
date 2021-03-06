# mri/z-test.jl

using Test: @testset

list = [
"exp_mult.jl"
"exp_xform.jl"
"coil_compress.jl"
"image_geom_mri.jl"
"kspace.jl"
"kspace-spiral.jl"
"mri_objects.jl"
"mri_trajectory.jl"
"sensemap-sim.jl"
]

for file in list
	@testset "$file" begin
		include(file)
	end
end
