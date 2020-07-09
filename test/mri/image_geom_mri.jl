# image_geom_mri.jl

using MIRT: image_geom_mri
using MIRT: MIRT_image_geom
using Test: @test, @inferred


@test (@inferred image_geom_mri( ; nx = 64, dx = 2)) isa MIRT_image_geom
