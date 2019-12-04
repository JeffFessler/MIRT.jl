# mri/z-test.jl

using Test: @test

@test ir_mri_coil_compress(:test)
@test ir_mri_kspace_ga_radial(:test)
@test ir_mri_sensemap_sim(:test)
@test image_geom_mri(:test)
@test mri_objects(:test)
