# image_geom_mri.jl

using MIRT: image_geom_mri
using MIRT: ImageGeom
using Test: @test, @inferred


#@inferred
@test image_geom_mri( ; nx = 64, dx = 2) isa ImageGeom{2}
