#=
rect_im.jl
2019-07, Helena H
=#

export rect_im

using Plots
# using MIRT: jim, MIRT_image_geom, downsample2


"""
`phantom = rect_im(ig, params;
	oversample=1, hu_scale=1, fov=maximum(ig.fovs), chat=false, how=:auto,
	replace=false, return_params=false)`

generate rectangle phantom image from parameters:
 `[x_center y_center x_width y_width angle_degrees amplitude]`

in
- `ig`				`image_geom()` object
- `params`			`[Nrect,6]` rect parameters. if empty use default

options
- `oversample::Integer`	oversampling factor, for grayscale boundaries
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
function rect_im(ig::MIRT_image_geom,
		params::AbstractArray{<:Real,2};
		oversample::Integer = 1,
		hu_scale::Real = 1,
		fov::Real = maximum(ig.fovs),
		chat::Bool = false,
		how::Symbol = :auto,
		replace::Bool = false,
		return_params::Bool = false)

	size(params,2) != 6 && throw("bad rect parameter vector size")

	if oversample > 1
		ig = ig.over(oversample)
	end
	args = (ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y, replace)

	params[:,6] .*= hu_scale

	do_fast = params[:,5] .== 0 # default for :auto
	if how == :fast
		do_fast[:] .= true
	elseif how == :slow
		do_fast[:] .= false
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
`phantom = rect_im_fast()`
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
		theta = rect[5] * (pi/180)
		theta != 0 && throw("theta=0 required")
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
`phantom = rect_im_slow()`
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
		theta = rect[5] * (pi/180) # Float64
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
`phantom = rect_im(nx, dx, params; args...)`

square image of size `nx` by `nx`,
specifying pixel size `dx` and rect `params`
"""
function rect_im(nx::Integer, dx::Real, params; args...)
	ig = image_geom(nx=nx, dx=1)
	return rect_im(ig, params; args...)
end


"""
`phantom = rect_im(nx::Integer, params; args...)`

square image of size `nx` by `nx` with
pixel size `dx=1` and rect `params`
"""
function rect_im(nx::Integer, params; args...)
	return rect_im(nx, 1., params; args...)
end


"""
`phantom = rect_im(nx::Integer; ny::Integer=nx, dx::Real=1)`

image of size `nx` by `ny` (default `nx`) with specified `dx` (default 1),
defaults to `:my_rect`
"""
function rect_im(nx::Integer; ny::Integer=nx, dx::Real=1,
	params::Symbol=:my_rect, args...)
	ig = image_geom(nx=nx, ny=ny, dx=dx)
	return rect_im(ig, params; args...)
end


"""
`phantom = rect_im(nx::Integer, ny::Integer; args...)`

`:my_rect` of size `nx` by `ny`
"""
function rect_im(nx::Integer, ny::Integer; args...)
	return rect_im(nx, ny=ny, dx=1.; args...)
end


"""
`phantom = rect_im(ig, code; args...)`

`code = :my_rect | :default | :smiley`
"""
function rect_im(ig::MIRT_image_geom, params::Symbol; args...)
	fov = ig.fovs
	if params == :my_rect
		params = my_rect(fov...)
	elseif params == :default
		params = rect_im_default_parameters(fov...)
	elseif params == :smiley
		params = smiley_parameters(fov...)
	else
		throw("bad phantom symbol $params")
	end
	return rect_im(ig, params; args...)
end


"""
`phantom = rect_im(ig; args...)`

`:default` (default) for given image geometry `ig`
"""
function rect_im(ig::MIRT_image_geom; args...)
	return rect_im(ig, :default; args...)
end



"""
`params = rect_im_default_parameters(xfov, yfov)`

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
`params = my_rect(xfov, yfov)`
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


"""
`rect_im()`

show docstring(s)
"""
function rect_im()
	@doc rect_im
end


"""
`rect_im_show()`
"""
function rect_im_show()
	#plot(1:10, 1:10)
	ig = image_geom(nx=2^8, ny=2^8, fov=100)

	x0 = rect_im(ig, [[[0.5, 0, 3, 20]*ig.dx..., 0, 1]';], oversample=3)
	p1 = jim(x0)

	x1 = rect_im(ig, :default; oversample=3, chat=true)
	p2 = jim(x1, title="default rects")

	x2 = rect_im(ig, :my_rect; oversample=3, chat=true)
	p3 = jim(x2, title="my rect")

	x3 =rect_im(ig, :smiley; oversample=3, chat=true)
	p4 = jim(x3, title="smiley")

	plot(p1, p2, p3, p4)
end


function rect_im_test()
	fov = 100
	rect_im_default_parameters(fov, fov)
	rect_im_show()
	true
end


"""
`rect_im(:test)`

`rect_im(:show)`

run tests
"""
function rect_im(test::Symbol)
	if test == :show
		return rect_im_show()
	end
	test != :test && throw(ArgumentError("test $test"))
	ig = image_geom(nx=2^8, dx=3)
	@test_throws String rect_im(ig, :bad)
	rect_im(ig, :smiley, how=:fast, replace=true)
	rect_im(ig, how=:slow, replace=true)
	rect_im(32, ny=30, dx=3, params=:default, how=:slow, return_params=true)
	rect_im(32, 30)
	@test_throws String rect_im(32, :default, how=:bad)
	rect_im(:show)
	@test rect_im_test()
	rect_im()
	true
end
