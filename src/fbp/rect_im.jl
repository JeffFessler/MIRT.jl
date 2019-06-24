using Plots

"""
`phantom = rect_im(ig, params;
oversample=1, hu_scale=1, fov=ig.fov, chat=false, how=:auto, replace=false, return_params=false)`

generate rectangle phantom image from parameters:
    `[x_center y_center x_width y_width angle_degrees amplitude]`

in
    `ig`			image_geom() object
    `params`		[Nrect,6] rect parameters

options
    `oversample`	int     oversampling factor, for grayscale boundaries
    `hu_scale`		float   use 1000 to scale
    `fov`           float   default ig.fov
	`chat`		   	bool
	`how`			symbol	:fast or :slow
	`replace`		bool
	`return_params` bool	if true, return both phantom and params

out
    `phantom`			[nx ny] image
	`params`			[Nrect,6] rect parameters (only return if return_params=true)

"""
function rect_im(ig::MIRT_image_geom,
    params::abstractArray{<:Real,2};
    oversample::Integer=1,
    hu_scale::Real=1,
    fov::Real=ig.fov,
	chat::Bool=false,
	how::Symbol=:auto,
	replace::Bool=false,
	return_params::Bool=false)

    args = (ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y, replace)

    if oversample > 1
    #    ig = ig.over(oversample)
		@warn("not yet done")
    end

	if isempty(fov)
		fov = ig.fov
	end

    params[:,6] .*= hu_scale

    do_fast = params[:,5] == 0 # default for :auto
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

	if any(!do_fast)
		phantom += rect_im_slow(params[!do_fast,:], args..., over)
	end

	if return_params
		return (phantom, params)
	end
	return phantom
end

"""
`phantom = rect_im_fast()`
"""
function rect_im_fast(params_in, nx, ny, dx, dy, offset_x, offset_y, replace)
    params = copy(params_in)
    if size(params,2) != 6
		throw("bad ellipse parameter vector size")
	end

	phantom = zeros(Float32, nx, ny)

	wx = (nx-1)/2 + offset_x
	wy = (ny-1)/2 + offset_y
	x1 = ([0:nx-1]' - wx) * dx
	y1 = ([0:ny-1] - wy) * dy

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
		if theta != 0
			throw("theta=0 required")
		end
		x = x1 - cx
		y = y1 - cy
		tx = fun(x-abs(dx)/2, x+abs(dx)/2, wx) / abs(dx)
		ty = fun(y-abs(dy)/2, y+abs(dy)/2, wy) / abs(dy)
		tmp = single(tx) * single(ty) # outer product (separable)
		if replace
			phantom(tmp > 0) = rect[6]
		else
			phantom = phantom + rect[6] * tmp
		end
	return phantom
end

"""
`phantom = rect_im_slow()`
"""
function rect_im_slow(params_in, ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y, replace, over)

	if size(params,2) != 6
		throw("bad rect parameter vector size")
	end

	phantom = zeros(Float32, nx*over, ny*over)

	wx = (nx*over - 1)/2 + offset_x*over
	wy = (ny*over - 1)/2 + offset_y*over
	xx = ([0:nx*over-1] - wx) * dx / over
	yy = ([0:ny*over-1] - wy) * dy / over
	(xx, yy) = ndgrid(xx, yy)

	# ticker reset
	ne = size(params)[1]
	for ie in 1:ne
		ticker(mfilename, ie, ne)

		rect = params[ie, :]
		cx = rect[1]
		wx = rect[3]
		cy = rect[2]
		wy = rect[4]
		theta = rect[5] * (pi/180)
		x = cos(theta) * (xx-cx) + sin(theta) * (yy-cy)
		y = -sin(theta) * (xx-cx) + cos(theta) * (yy-cy)
		tmp = abs(x / wx) < 1/2 & abs(y / wy) < 1/2 # ?
		if replace
			phantom[tmp > 0] = rect[6] # ?
		else
			phantom = phantom + rect[6] * tmp
		end
	end
	return phantom
end

function rect_im_default_parameters(xfov, yfov)
	f = 1/64
	params = [
		0     0     50     50     0     1
		10    -16	25	   16	  0	    -0.5
		-13	  15	13	   13	  1*45  1
		-18   0     1      1      0     1
		-12   0     1      1      0     1
		-6    0     1      1      0     1
		0     0     1      1      0     1
		6     0     1      1      0     1
		12    0     1      1      0     1
		18    0     1      1      0     1
	]

	params[:,[1,3]] .*= xfov/64
	params[:,[2,4]] .*= yfov/64

	return params
end



#=
if isempty(params)
	params = rect_im_default_parameters(xfov, yfov)
end
=#

"""
`(xx,yy) = ndgrid(x,y)`
"""
function ndgrid(x::AbstractVector{<:Number},
				y::AbstractVector{<:Number})
	return (repeat(x, 1, length(y)), repeat(y', length(x), 1))
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
function rect_im_show() # ?
	#plot(1:10, 1:10)
	ig = image_geom(nx=2^8, ny=2^8, fov=100) # ?
	x0 = rect_im(ig, [[0.5, 0, 3, 20].*ig.dx, 0, 1], oversample=3)
	p1 = jim(x0)
	x0 = rect_im(ig, [], oversample=3, chat=true, im)
	p2 = jim(ig.x, ig.y, x0, title="default rects")

	plot(p1, p2)
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
	rect_im()
	rect_im(:show)
	rect_im_test()
end
