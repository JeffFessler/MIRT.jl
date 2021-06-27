#=
rect_im.jl
2019-07, Helena H
=#

export rect_im

#using MIRT: ImageGeom, downsample2


"""
    phantom = rect_im(ig, params ;
	oversample=1, hu_scale=1, fov=maximum(ig.fovs), chat=false, how=:auto,
	replace=false, return_params=false)

generate rectangle phantom image from parameters:
 `[x_center y_center x_width y_width angle_degrees amplitude]`

in
- `ig`				`image_geom()` object
- `params`			`[Nrect,6]` rect parameters. if empty use default

options
- `oversample::Int`	oversampling factor, for grayscale boundaries
- `hu_scale::Real`		use 1000 to scale
- `fov::Real`			default `maximum(ig.fovs)`
- `chat::Bool`			verbosity?
- `how::Symbol`			`:fast` or `:slow`; default `:auto`
- `replace::Bool`		default `false`
- `return_params::Bool`	if true, return both phantom and params

out
- `phantom`		`[nx ny]` image (Float32)
- `params`		`[Nrect 6]` rect parameters (only return if `return_params=true`)
"""
function rect_im(
	ig::ImageGeom,
	params::AbstractArray{<:Real,2} ;
	oversample::Int = 1,
	hu_scale::Real = 1,
	fov::Real = maximum(ig.fovs),
	chat::Bool = false,
	how::Symbol = :auto,
	replace::Bool = false,
	return_params::Bool = false,
)

	size(params,2) != 6 && throw("bad rect parameter vector size")

	if oversample > 1
		ig = ig.over(oversample)
	end
	args = (ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y, replace)

	params[:,6] .*= hu_scale

	do_fast = params[:,5] .== 0 # default for :auto
	if how === :fast
		vec(do_fast) .= true
	elseif how === :slow
		vec(do_fast) .= false
	elseif how != :auto
		throw("bad how :how")
	end

	phantom = zeros(Float32, ig.nx, ig.ny)

	if any(do_fast)
		phantom += rect_im_fast(params[do_fast,:], args...)
	end

	if any(.!do_fast)
		tmp = rect_im_slow(params[.!do_fast,:], args...)
		phantom += tmp
	end

	if oversample > 1
		phantom = downsample2(phantom, oversample)
	end

	if return_params
		return (phantom, params)
	end

	return phantom
end


"""
phantom = rect_im_fast()
for non-rotated rectangles
	using exact integration over each pixel so over-sampling is irrelevant
"""
function rect_im_fast(params, nx, ny, dx, dy, offset_x, offset_y, replace)
	phantom = zeros(Float32, nx, ny)

	wx = (nx-1)/2 + offset_x # scalars
	wy = (ny-1)/2 + offset_y
	x1 = ((0:nx-1) .- wx) * dx # arrays
	y1 = ((0:ny-1) .- wy) * dy
	fun = (x1, x2, wx) -> # integrated rect(x/wx) function from x1 to x2
		max(min(x2, wx/2) - max(x1, -wx/2), 0)

	# ticker reset
	ne = size(params)[1]
	for ie in 1:ne
		#ticker(mfilename, ie, ne)
		rect = params[ie, :]
		cx = rect[1]
		wx = rect[3]
		cy = rect[2]
		wy = rect[4]
		theta = deg2rad(rect[5])
		theta != zero(theta) && throw("theta=0 required")
		value = Float32(rect[6])

		x = x1 .- cx
		y = y1 .- cy
		tx = fun.(x .- abs(dx)/2, x .+ abs(dx)/2, wx) / abs(dx)
		ty = fun.(y .- abs(dy)/2, y .+ abs(dy)/2, wy) / abs(dy)
		tmp = Float32.(tx) * Float32.(ty)' # outer product (separable)
		if replace
			phantom[tmp .> 0] .= value
		else
			phantom += value * tmp
		end
	end
	return phantom
end


"""
    phantom = rect_im_slow()
for rotated rectangles
"""
function rect_im_slow(params_in, nx, ny, dx, dy, offset_x, offset_y, replace)

	params = copy(params_in)
	phantom = zeros(Float32, nx, ny)

	wx = (nx - 1)/2 + offset_x # scalar
	wy = (ny - 1)/2 + offset_y
	x = ((0:nx-1) .- wx) * dx
	y = ((0:ny-1) .- wy) * dy
	(xx, yy) = ndgrid(x, y)

	# ticker reset
	ne = size(params)[1]
	for ie in 1:ne
		#ticker(mfilename, ie, ne)

		rect = params[ie, :]
		cx = rect[1]
		wx = rect[3]
		cy = rect[2]
		wy = rect[4]
		theta = deg2rad(rect[5])
		value = Float32(rect[6])

		(x,y) = rotate2d(xx .- cx, yy .- cy, theta)
		tmp = (abs.(x / wx) .< 1/2) .& (abs.(y / wy) .< 1/2) # 2d rect

		if replace
			phantom[tmp .> 0] .= value
		else
			phantom += value * tmp
		end
	end
	return phantom
end


"""
    phantom = rect_im(nx, dx, params ; args...)

square image of size `nx` by `nx`,
specifying pixel size `dx` and rect `params`
"""
function rect_im(nx::Int, dx::Real, params ; args...)
	ig = image_geom(nx=nx, dx=1)
	return rect_im(ig, params ; args...)
end


"""
    phantom = rect_im(nx::Int, params ; args...)

square image of size `nx` by `nx` with
pixel size `dx=1` and rect `params`
"""
function rect_im(nx::Int, params ; args...)
	return rect_im(nx, 1., params ; args...)
end


"""
    phantom = rect_im(nx::Int ; ny::Int=nx, dx::Real=1)

image of size `nx` by `ny` (default `nx`) with specified `dx` (default 1),
defaults to `:my_rect`
"""
function rect_im(nx::Int ; ny::Int=nx, dx::Real=1,
	params::Symbol=:my_rect, args...)
	ig = image_geom(nx=nx, ny=ny, dx=dx)
	return rect_im(ig, params ; args...)
end


"""
    phantom = rect_im(nx::Int, ny::Int ; args...)

`:my_rect` of size `nx` by `ny`
"""
function rect_im(nx::Int, ny::Int ; args...)
	return rect_im(nx, ny=ny, dx=1. ; args...)
end


"""
    phantom = rect_im(ig, code ; args...)

`code = :my_rect | :default | :smiley`
"""
function rect_im(ig::ImageGeom, params::Symbol ; args...)
	fov = ig.fovs
	if params === :my_rect
		params = my_rect(fov...)
	elseif params === :default
		params = rect_im_default_parameters(fov...)
	elseif params === :smiley
		params = smiley_parameters(fov...)
	else
		throw("bad phantom symbol $params")
	end
	return rect_im(ig, params ; args...)
end


"""
    phantom = rect_im(ig ; args...)

`:default` (default) for given image geometry `ig`
"""
rect_im(ig::ImageGeom ; args...) = rect_im(ig, :default ; args...)



"""
params = rect_im_default_parameters(xfov, yfov)

default parameters
"""
function rect_im_default_parameters(xfov, yfov)
	f = 1/64
	params = [
		0	0	50	50	0	1
		10	-16	25	16	0	-0.5
		-13	15	13	13	1*45	1
		-18	0	1	1	0	1
		-12	0	1	1	0	1
		-6	0	1	1	0	1
		0	0	1	1	0	1
		6	0	1	1	0	1
		12	0	1	1	0	1
		18	0	1	1	0	1
	]

	params[:,[1,3]] .*= xfov/64 # x_center and x_width
	params[:,[2,4]] .*= yfov/64 # y_center and y_width

	return params
end


"""
params = my_rect(xfov, yfov)
"""
function my_rect(xfov, yfov)
	f = 1/64
	rect = [
		35	35	20	10	45	1
		35	-35	20	10	-45	1
		-35	-35	20	10	45	1
		-35	35	20	10	-45	1
		0	0	40	40	45	1
		0	0	40	40	0	0.5
		40	0	12	1.5	0	1.5
		-40	0	12	1.5	0	1.5
	]
	return rect
end


"""
`params = smiley_parameters(xfov, yfov)`

smiley face out of rects
"""
function smiley_parameters(xfov, yfov)
	rect = [
		0	0	80	80	0	0.5
		-20	-20	10	15	0	1 #eyes
		20	-20	10	15	0	1
		0	25	50	6	0	1 #mouth
		-23	19	4	6	0	1
		23	19	4	6	0	1
	]
	return rect
end
