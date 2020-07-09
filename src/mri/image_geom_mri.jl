#=
Project:      MIRTdev
File:         image_geom_mri.jl
Description:  Image geom function
Author:       Rajas Gupta
Date Created: October 1st 2019
Date Update:  October 7th 2019

2007-02-22, Jeff Fessler, University of Michigan (image_geom_mri.m)
=#

export image_geom_mri

#using MIRT: image_geom, MIRT_image_geom


"""
    ig = image_geom_mri(varargin)

Same as `image_geom()` but default offsets are 0.5 so that
image pixel indices go from -N/2 to N/2-1.

"""
function image_geom_mri( ; kwargs...)
    return image_geom( ; offset_x = .5, offset_y = .5, offset_z = .5, kwargs...)
end
