# mri/z-list.jl

include("exp_mult.jl")
include("exp_xform.jl")
include("coil_compress.jl")
include("kspace.jl")
include("kspace-spiral.jl")
include("mri_objects.jl")
include("mri_trajectory.jl")
include("sensemap-sim.jl")

function image_geom_mri(args...; kwargs...)
    throw("image_geom_mri is deprecated; use ImageGeoms with 'offsets=:dsp'")
end
