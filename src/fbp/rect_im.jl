using Plots

"""
generate rectangle phantom image from parameters:
    [x_center y_center x_width y_width angle_degrees amplitude]

in
    ig          image_geom() object
    params      [Nrect,6] rect parameters

options
    `oversample'    int     oversampling factor, for grayscale boundaries
    'hu_scale'      float   use 1000 to scale
    'fov'           float   default ig.fov
	chat
	how

out
    phantom [nx ny] image

"""
function rect_im(ig::MIRT_image_geom,
    params::abstractArry{<:Real,2};
    oversample::Integer=1,
    hu_scale::Real=1,
    fov::Real=ig.fov,
	chat::Bool=false,
	how::Symbol=:auto,
	)

    args = (ig, params, oversample, hu_scale, fov, chat, how)

    if oversample > 1
    #    ig = ig.over(oversample)
		@warn("not yet done")
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
    
    phantom = zeros(Float32, ig.nx, ig.,ny)
    
    if any(do_fast)
        phantom += rect_im_fast(params todo)
	end

	if any(!.do_fast)
		phantom += rect_im_slow(params)
	end

	return phantom
end


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

	# not done
	#fun = @(x1, x2, wx) ... % integrated rect(x/wx) function from x1 to x2
	#max(min(x2, wx/2) - max(x1, -wx/2), 0);

	#ticker reset
	ne = size(params)[1]
	for ie in 1:ne
		#ticker(mfilename, ie, ne)

		rect = params[ie, :]
		cx = rect[1]
		wx = rect[3]
		cy = rect[2]
		wy = rect[4]
		theta = rect[5]
		if theta ~= 0, fail 'theta=0 required', end
		x = x1 - cx;
		y = y1 - cy;
		tx = fun(x-abs(dx)/2, x+abs(dx)/2, wx) / abs(dx);
		ty = fun(y-abs(dy)/2, y+abs(dy)/2, wy) / abs(dy);
		tmp = single(tx) * single(ty); % outer product (separable)
		if replace
			phantom(tmp > 0) = rect(6);
		else
			phantom = phantom + rect(6) * tmp;
		end
end
#=



function rect_im_slow(params_in, nx, ny, dx, dy, offset_x, offset_y, over, replace)

	if size(params,2) != 6
		throw("bad rect parameter vector size")
	end

	phantom = zeros(Float32, nx*over, ny*over)

	wx = (nx*over - 1)/2 + offset_x*over
	wy = (ny*over - 1)/2 + offset_y*over
	xx = ([0:nx*over-1] - wx) * dx / over
	yy = ([0:ny*over-1] - wy) * dy / over
	# [xx, yy] = ndgrid(xx, yy)

end =#



#=
if isempty(params)
	params = rect_im_default_parameters(xfov, yfov)
end
=#
