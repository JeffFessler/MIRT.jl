# disk-phantom.jl

using MIRT: disk_phantom_params, ellipse_im, jim
using Test: @inferred


ig = image_geom(nx = 128, fov = 240)
params = @inferred disk_phantom_params( ; minsep = ig.dx*8)
@test params isa Matrix{Float32}

tmp = ellipse_im(ig, params, oversample=3)
jim(ig.x, ig.y, tmp)
