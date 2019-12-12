# mri/z-test.jl

using Test: @test

@test exp_xform(:test)
@test image_geom_mri(:test)
@test ir_mri_coil_compress(:test)
@test ir_mri_kspace_ga_radial(:test)
@test ir_mri_sensemap_sim(:test)
@test mri_kspace_spiral(:test)
@test mri_objects(:test)
