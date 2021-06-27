#=
ellipse_im.jl
Copyright 2019-03-05, Jeff Fessler, University of Michigan
=#

export ellipse_im, ellipse_im_params

# using MIRT: image_geom, disk_phantom_params, downsample2, rotate2d


"""
    phantom = ellipse_im(ig, params ;
	rot=0, oversample=1, hu_scale=1, replace=false)

Generate ellipse phantom image from parameters:

`[x_center y_center x_radius y_radius angle_degrees amplitude]`

in
- `ig`		from `image_geom()`
- `params`	`[ne 6]` ellipse parameters

# Arguments
- `rot`			rotate ellipses by this amount [degrees]
- `oversample`	oversampling factor, for grayscale boundaries
- `hu_scale`	use 1000 to scale shepp-logan to HU
- `replace`		replace ellipse values if true, else add
- `how`			`:fast` is the only option

out
- `phantom`		[nx ny]	image (Float32)

note: `op ellipse` in aspire with `nsub=3` is `oversample=4 = 2^(3-1)` here

"""
function ellipse_im(
	ig::ImageGeom,
	params::AbstractMatrix{<:RealU} ;
	rot::RealU = 0,
	oversample::Int = 1,
	hu_scale::RealU = 1,
	replace::Bool = false,
	how::Symbol = :fast, # todo
)

    size(params,2) == 6 || throw("bad ellipse parameter vector size")

#=
	if how === :fast && oversample == 1
		warning("ignoring :fast option for oversample=1")
		how = :slow
	end
=#

	if oversample > 1
		ig = ig.over(oversample)
	end

	args = (ig.nx, ig.ny, params, ig.dx, ig.dy, ig.offset_x, ig.offset_y,
        rot, oversample, replace, hu_scale)

    T = promote_type(Float32, typeof.(params[:,6])...) # in case of units
    phantom = zeros(T, ig.dims)

	if how === :fast
        ellipse_im_fast!(phantom, args...)
#	elseif how === :slow
#       ellipse_im_slow!(phantom, args...)
	else
		throw("bad how $how")
	end

	if oversample > 1
		phantom = downsample2(phantom, oversample)
	end

	return phantom
end


#=
%
% ellipse_im_slow()
% brute force fine grid - can use lots of memory
%
function phantom = ellipse_im_slow(nx, ny, params, dx, dy, ...
	#offset_x, offset_y, rot, over, replace)

if size(params,2) ~= 6
	fail('bad ellipse parameter vector size')
end

% optional rotation of ellipse parameters
if rot ~= 0
	th = deg2rad(rot)
	cx = params[:,1)
	cy = params[:,2)
	params[:,1) = cx * cos(th) + cy * sin(th)
	params[:,2) = -cx * sin(th) + cy * cos(th)
	params[:,5) = params[:,5) + rot
	clear cx cy th
end

wx = (nx*over-1)/2 + offset_x * over
wy = (ny*over-1)/2 + offset_y * over
xx = ((0:nx*over-1) - wx) / over * dx
yy = ((0:ny*over-1) - wy) / over * dy
[xx yy] = ndgrid(xx, yy); % fine grid, equally spaced

phantom = zeros(nx*over, ny*over, 'single'); % fine array

ticker reset
ne = nrow(params)
for ie = 1:ne
	ticker(mfilename, ie, ne)

	ell = params[ie, :)
	cx = ell(1);	rx = ell(3)
	cy = ell(2);	ry = ell(4)
	theta = deg2rad(ell(5))
	[xr yr] = rot2(xx-cx, yy-cy, theta)
	tmp = (xr / rx).^2 + (yr / ry).^2 <= 1

	if replace
		phantom(tmp > 0) = ell(6)
	else
		phantom = phantom + ell(6) * tmp
	end
end

phantom = downsample2(phantom, over)
=#


"""
phantom = ellipse_im_fast!()
"""
function ellipse_im_fast!(phantom, nx, ny, params_in, dx, dy,
    offset_x, offset_y, rot, over, replace, hu_scale,
)

    params = copy(params_in)

	# optional rotation
	if rot != zero(rot)
        th = deg2rad(rot)
		cx = params[:,1]
		cy = params[:,2]
		params[:,1] = cx * cos(th) + cy * sin(th)
		params[:,2] = -cx * sin(th) + cy * cos(th)
		params[:,5] .+= rot
	end

	wx = (nx-1)/2 + offset_x
	wy = (ny-1)/2 + offset_y
	x1 = ((0:nx-1) .- wx) * dx
	y1 = ((0:ny-1) .- wy) * dy
	(xx, yy) = ndgrid(x1, y1)

	hx = abs(dx) / 2
	hy = abs(dy) / 2

	for ie = 1:size(params,1)

		ell = params[ie, :]
		cx = ell[1]; rx = ell[3]
		cy = ell[2]; ry = ell[4]
		theta = deg2rad(ell[5])
        value = ell[6] * hu_scale

		xs = xx .- cx # shift per ellipse center
		ys = yy .- cy

		# coordinates of "outer" corner of each pixel, relative to ellipse center
		xo = xs + sign.(xs) * hx
		yo = ys + sign.(ys) * hy

		# voxels that are entirely inside the ellipse:
		(xr, yr) = rotate2d(xo, yo, theta)
		is_inside = (xr / rx).^2 + (yr / ry).^2 .<= 1

		if replace
			phantom[is_inside] .= value
		else
            phantom .+= value * is_inside
		end

	end # ie loop

	return phantom
end # ellipse_im_fast()


"""
    phantom = ellipse_im(nx, dx, params ; kwarg...)

square image of size `nx × nx`,
specifying pixel size `dx` and ellipse `params`
"""
function ellipse_im(nx::Int, dx::RealU, params ; kwarg...)
    ig = image_geom( ; nx, dx)
	return ellipse_im(ig, params ; kwarg...)
end


"""
    phantom = ellipse_im(nx::Int, params ; kwarg...)

square image of size `nx` by `nx` with
pixel size `dx=1` and ellipse `params`
"""
function ellipse_im(nx::Int, params ; kwarg...)
	return ellipse_im(nx, 1., params ; kwarg...)
end


"""
    phantom = ellipse_im(nx::Int ; ny::Int=nx, dx::Number=1, kwarg...)

image of size `nx` by `ny` (default `nx`) with specified `dx` (default 1),
defaults to `:shepplogan_emis`
"""
function ellipse_im(nx::Int ; ny::Int=nx, dx::RealU=1, kwarg...)
    ig = image_geom( ; nx, ny, dx)
	return ellipse_im(ig, :shepplogan_emis ; kwarg...)
end


"""
    phantom = ellipse_im(nx::Int, ny::Int ; kwarg...)

`:shepplogan_emis` of size `nx` by `ny`
"""
function ellipse_im(nx::Int, ny::Int ; kwarg...)
    return ellipse_im(nx ; ny, dx=1., kwarg...)
end


"""
    phantom = ellipse_im(ig, code ; kwarg...)

`code = :shepplogan | :shepplogan_emis | :shepplogan_brainweb | :southpark`
"""
function ellipse_im(ig::ImageGeom, params::Symbol ; kwarg...)
	params = ellipse_im_params(ig, params)
	return ellipse_im(ig, params ; kwarg...)
end


"""
    phantom = ellipse_im(ig ; kwarg...)

`:shepplogan` (default) for given image geometry `ig`
"""
function ellipse_im(ig::ImageGeom ; kwarg...)
	return ellipse_im(ig, :shepplogan; kwarg...)
end


"""
    params = ellipse_im_params(ig::ImageGeom, params::Symbol)

`code = :shepplogan | :shepplogan_emis | :shepplogan_brainweb | :southpark | :disks`
"""
function ellipse_im_params(ig::ImageGeom, params::Symbol)
	if params === :disks
		params = disk_phantom_params(fov=ig.fovs[1])
	elseif params === :shepplogan || params === :kak
		params = shepp_logan_parameters(ig.fovs..., case=:kak)
	elseif params === :shepplogan_emis || params === :emis
		params = shepp_logan_parameters(ig.fovs..., case=:emis)
	elseif params === :shepplogan_brainweb || params === :brainweb
		params = shepp_logan_parameters(ig.fovs..., case=:brainweb)
	elseif params === :southpark
		params = south_park_parameters(fov=ig.fovs[1])
	else
		throw("bad phantom symbol $params")
	end
	return params
end


"""
params = shepp_logan_parameters(xfov, yfov ; case::Symbol ; u::Tuple)

Shepp Logan phantom parameters from Kak and Slaney text, p. 255.

By default the first four columns are unitless "fractions of field of view",
so columns 1,3 are scaled by `xfov` and columns 2,4 are scaled by `xfov`.
The optional 3-tuple `u` specifies scaling and/or units:
* columns 1-4 (center, radii) are scaled by `u[1]` (e.g., mm),
* column 5 (angle) is scaled by `u[2]` (e.g., `1` or `°`),
* column 6 (value) is scaled by `u[3]` (e.g., `1/cm`) for an attenuation map.
"""
function shepp_logan_parameters(xfov::RealU, yfov::RealU ;
    case::Symbol=:kak,
    u::Tuple{a,b,c} where {a,b,c} = (1,1,1), # unit scaling
)
	params = [
	0		0		0.92	0.69	90	2
	0		-0.0184	0.874	0.6624	90	-0.98
	0.22	0		0.31	0.11	72	-0.02
	-0.22	0		0.41	0.16	108	-0.02
	0		0.35	0.25	0.21	90	0.01
	0		0.1		0.046	0.046	0	0.01
	0		-0.1	0.046	0.046	0	0.01
	-0.08	-0.605	0.046	0.023	0	0.01
	0		-0.605	0.023	0.023	0	0.01
	0.06	-0.605	0.046	0.023	90	0.01]

	params[:,[1,3]] .*= xfov/2
	params[:,[2,4]] .*= yfov/2

	if case === :emis
		params[:,6] = [1, 1, -2, 2, 3, 4, 5, 6, 1, 1]
	elseif case === :brainweb
		params[:,6] = [1, 0, 2, 3, 4, 5, 6, 7, 8, 9] # brainweb uses index 1-10
	elseif case != :kak
		throw("bad phantom case $case")
	end

    return [params[:,1:4] * u[1] params[:,5] * u[2] params[:,6] * u[3]]
end


"""
param = south_park_parameters( ; fov::Number = 100)
"""
function south_park_parameters( ; fov::RealU = 100)
	xell = [
		0. 0 85 115 0 100
		0 -60 30 20 0 -80 # mouth
		30 20 25 35 30 20 # eyes
		-30 20 25 35 -30 20
		35 25 7 7 0 -100 # pupils
		-15 25 7 7 0 -100
		0 75 60 15 0 -50] # hat
	xell[:,1:4] .*= fov/256
	return xell
end
