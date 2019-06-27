# mri/z-test.jl

using Test

@test ir_mri_kspace_ga_radial(:test)
@test ir_mri_sensemap_sim(:test)
