# fbp/z-list.jl

include("sino_geom.jl")
include("sino_plot.jl") # todo: move to democard

include("rotate2d.jl")

include("cuboid_im.jl") # must be after image_geom.jl
include("ellipsoid_im.jl")
include("rect_im.jl")


# deprecations below here:

export image_geom
function image_geom(args...; kwargs...)
    throw("image_geom is deprecated; use ImageGeoms")
end

export cuboid_im
function cuboid_im(args...; kwargs...)
    throw("cuboid_im is deprecated; use ImagePhantoms")
end


export disk_phantom_params
function disk_phantom_params(args...; kwargs...)
    throw("disk_phantom_params is deprecated; use ImagePhantoms")
end


export ellipse_im, ellipse_im_params
function ellipse_im(args...; kwargs...)
    throw("ellipse_im is deprecated; use ImagePhantoms")
end
function ellipse_im_params(args...; kwargs...)
    throw("ellipse_im_params is deprecated; use ImagePhantoms")
end

export ellipse_sino
function ellipse_sino(args...; kwargs...)
    throw("ellipse_sino is deprecated; use ImagePhantoms")
end

export rect_sino
function rect_sino(args...; kwargs...)
    throw("rect_sino is deprecated; use ImagePhantoms")
end
