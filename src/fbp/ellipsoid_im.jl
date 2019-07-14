#=
ellipsoid_im.jl
2019-07 Helena H
=#

export ellipsoid_im

using Plots
using Test: @test, @test_throws
#using MIRT: jim, MIRT_image_geom, image_geom


"""
`phantom = ellipsoid_im(ig, params;
	oversample=1, checkfov=false, how=:slow, showmem=false, hu_scale=1, return_params=false)`

generate ellipsoid phantom image from parameters:
	`[x_center, y_center, z_center, x_radius, y_radius, z_radius,
		xy_angle_degrees, z_angle_degrees, density]`

in
- `ig`			from `image_geom()`
- `params`		`[N 9]` ellipsoid parameters.

option
- `oversample::Integer`		oversampling factor (default:1)
- `checkfov::Bool`			warn if any ellipsoid is out of fov
- `how::Symbol`				`:fast` does it fast -- to do, only works slow
							`:lowmem` uses less memory than :fast but slower
							`:slow` default
- `showmem::Bool`
- `hu_scale::Real`			use 1000 to scale shepp-logan to HU; default 1
- `return_params::Bool`		if true, return both phantom and params

out
- `phantom`		`[nx ny nz]` image
- `params`		`[N 9]` ellipsoid parameters (only if `return_params=true`)
"""
function ellipsoid_im(ig::MIRT_image_geom, params::AbstractMatrix{<:Real} ;
		oversample::Integer = 1,
		checkfov::Bool = false,
		how::Symbol = :slow,
		showmem::Bool = false,
		hu_scale::Real = 1,
		return_params::Bool = false)

	size(params,2) != 9 && throw("bad cuboid parameter vector size")
	params[:,9] .*= hu_scale

	if oversample > 1
		ig = ig.over(oversample)
	end

	args = (ig.nx, ig.ny, ig.nz, params, ig.dx, ig.dy, ig.dz,
			ig.offset_x, ig.offset_y, ig.offset_z)
	phantom = zeros(Float32, ig.nx, ig.ny, ig.nz)

	checkfov && !ellipsoid_im_check_fov(args...) && throw("ellipsoid exceeds FOV")

	if how == :slow
		phantom += ellipsoid_im_slow(args..., showmem)
#=
	elseif how == :lowmem
		phantom += ellipsoid_im_lowmem(args..., showmem)
	elseif how == :fast
		phantom += ellipsoid_im_fast(args..., showmem)
=#
	else
		throw("bad how $how")
	end

	if oversample > 1
		phantom = downsample3(phantom, oversample)
	end

	if return_params
		return (phantom, params)
	end
	return phantom
end


"""
`ellipsoid_im_slow()`

brute force fine grid - can use lots of memory
"""
function ellipsoid_im_slow(nx, ny, nz, params, dx, dy, dz,
		offset_x, offset_y, offset_z, showmem)

	phantom = zeros(Float32, nx, ny, nz)

	wx = (nx - 1)/2 + offset_x
	wy = (ny - 1)/2 + offset_y
	wz = (nz - 1)/2 + offset_z
	x = ((0:(nx-1)) .- wx) * dx
	y = ((0:(ny-1)) .- wy) * dy
	z = ((0:(nz-1)) .- wz) * dz
	xmax = maximum(x)
	xmin = minimum(x)
	ymax = maximum(y)
	ymin = minimum(y)
	zmax = maximum(z)
	zmin = minimum(z)
	(xx, yy, zz) = ndgrid(x, y, z)

	#ticker reset
	np = size(params)[1]
	for ip in 1:np
		#ticker(ip, np)

		par = params[ip, :]
		cx = par[1]
		cy = par[2]
		cz = par[3]
		rx = par[4]
		ry = par[5]
		rz = par[6]

		azim = deg2rad(par[7])
		polar = deg2rad(par[8])
		(xr, yr, zr) = rot3(xx .- cx, yy .- cy, zz .- cz, azim, polar)
		tmp = ((xr / rx).^2 + (yr / ry).^2 + (zr / rz).^2) .<= 1
		value = Float32(par[9])
		phantom += value * tmp
	end
	return phantom
end


#=
"""
`ellipsoid_im_fast()`

currently not working
only slow option works
"""
function ellipsoid_im_fast(nx, ny, nz, params, dx, dy, dz,
	offset_x, offset_y, offset_z, over, showmem)

	phantom = zeros(Float32, nx, ny, nz)

	wx = (nx-1)/2 + offset_x
	wy = (ny-1)/2 + offset_y
	wz = (nz-1)/2 + offset_z
	x = ((0:nx-1) .- wx) * dx
	y = ((0:ny-1) .- wy) * dy
	z = ((0:nz-1) .- wz) * dz
	xmax = maximum(x)
	xmin = minimum(x)
	ymax = maximum(y)
	ymin = minimum(y)
	zmax = maximum(z)
	zmin = minimum(z)
	(xx, yy, zz) = ndgrid(x, y, z)

	if over > 1
		tmp = ((1:over) - (over+1)/2) / over
		(xf, yf, zf) = ndgrid(tmp*dx, tmp*dy, tmp*dz)
		xf = xf[:]'
		yf = yf[:]'
		zf = zf[:]'
	end

	hx = abs(dx) / 2
	hy = abs(dy) / 2
	hz = abs(dz) / 2

	#ticker reset
	np = size(params)[1]
	for ip in 1:np
		#ticker(mfilename, ip, np)
		par = params[ip,:]
		cx = par[1]
		cy = par[2]
		cz = par[3]
		rx = par[4]
		ry = par[5]
		rz = par[6]

		azim = par[7] * (pi/180)
		polar = par[8] * (pi/180)

		xs = xx .- cx
		ys = yy .- cy
		zs = zz .- cz

		(xr, yr, zr) = rot3(xs, ys, zs, azim, polar)
		if over == 1
			vi = ((xr / rx).^2 + (yr / ry).^2 + (zr / rz).^2) .<= 1
			value = Float32(par[9])
			phantom += value * vi
		end

		#= check all 8 corners of the cube that bounds
		all possible 3D rotations of the voxel =#
		hh = sqrt(hx^2 + hy^2 + hz^2)
		vi = true # voxels entirely inside the ellipsoid
		#vo = true # voxels entirely outside the ellipsoid

		for ix in -1:1, iy in -1:1, iz in -1:1
			xo = xr .+ ix*hh
			yo = yr .+ iy*hh
			zo = zr .+ iz*hh
			vi = vi .& (((xo/rx).^2 + (yo/ry).^2 + (zo/rz).^2) .< 1)
		end

		vo = false #to do. for now, "outside" test is failing
		if ip == 1
			throw("todo:must debug this")
		end
		if any(vi[:] & vo[:])
			throw("bug")
		end
		#=
		if 0
			% coordinates of "outer" corner of each voxel, relative to ellipsoid center
			sign_p = @(x) (x >= 0) * sqrt(3); % modified sign()
			xo = xr + sign_p(xr) * hx;
			yo = yr + sign_p(yr) * hy;
			zo = zr + sign_p(zr) * hz;

			% voxels that are entirely inside the ellipsoid
			vi = (xo / rx).^2 + (yo / ry).^2 + (zo / rz).^2 <= 1;
		end

		if 0
			% coordinates of "inner" corner of each pixel, relative to ellipse center
			sign_n = @(x) (x > 0) * sqrt(3); % modified sign()
			xi = xr - sign_n(xs) * hx;
			yi = yr - sign_n(ys) * hy;
			zi = zr - sign_n(zs) * hz;

			% voxels that are entirely outside the ellipsoid
			vo = (max(abs(xi),0) / rx).^2 + (max(abs(yi),0) / ry).^2 ...
				+ (max(abs(zi),0) / rz).^2 >= 1;
		end
		=#

		edge = !vi & !vo
		if showmem
			println("edge fraction ",
				sum(edge[:]) / prod(size(edge)), "=",
				sum(edge[:]), "/", prod(size(edge)))
		end

		x = xx[edge] .- cx
		y = yy[edge] .- cy
		z = zz[edge] .- cz
		x = x .+ xf'
		y = y .+ yf'
		z = z .+ zf'
		(xr, yr, zr) = rot3(x, y, z, azim, polar)
		in = ((xr/rx).^2 + (yr/ry).^2 + (zr/rz).^2) .<= 1
		tmp = mean(in, 2)
		gray = Float32(vi)
		gray[edge] = tmp
		#=
		if 0 % todo: help debug
			ee = (gray > 0) & (gray < 1);
			im(ee)
			im(vi)
		prompt
			%minmax(tmp)
			%clf, im pl 1 2
			%im(tmp)
		end=#

		value = Float32(par[9])
		phantom += value*gray
	end
	return phantom
end


"""
`ellipsoid_im_lowmem()`

Do 'one slice at a time' to reduce memory
"""
function ellipsoid_im_lowmem(nx, ny, nz, params, dx, dy, dz,
		offset_x, offset_y, offset_z)

	phantom = zeros(Float32, nx, ny, nz)
	for iz in 1:nz
		offset_z_new = (nz-1)/2 + offset_z - (iz-1)
		phantom[:,:,iz] = ellipsoid_im_fast(nx, ny, 1, params, dx, dy, dz,
		offset_x, offset_y, offset_z_new, showmem && iz == 1)
	end
	return phantom
end
=#


"""
`ellipsoid_im_check_fov()`
"""
function ellipsoid_im_check_fov(nx, ny, nz, params,
		dx, dy, dz, offset_x, offset_y, offset_z)
	wx = (nx - 1)/2 + offset_x
	wy = (ny - 1)/2 + offset_y
	wz = (nz - 1)/2 + offset_z
	xx = ((0:nx - 1) .- wx) * dx
	yy = ((0:ny - 1) .- wy) * dy
	zz = ((0:nz - 1) .- wz) * dz

	xmax = maximum(xx)
	xmin = minimum(xx)
	ymax = maximum(yy)
	ymin = minimum(yy)
	zmax = maximum(zz)
	zmin = minimum(zz)

	for ip in 1:size(params)[1]
		par = params[ip, :]
		cx = par[1]
		cy = par[2]
		cz = par[3]
		rx = par[4]
		ry = par[5]
		rz = par[6]

		if (cx + rx > xmax) || (cx - rx < xmin)
			@warn("fov: x range $xmin $xmax, cx=$cx, rx=$rx")
			return false
		end
		if (cy + ry > ymax) || (cy - ry < ymin)
			@warn("fov: y range $ymin $ymax, cy=$cy, ry=$ry")
			return false
		end

		if (cz + rz > zmax) || (cz - rz < zmin)
			@warn("fov: z range $zmin $zmax, cz=$cz, rz=$rz")
			return false
		end
	end
	return true
end


"""
`phantom = ellipsoid_im(ig, ptype ; args...)`

`ptype = :zhu | :kak | :e3d | :spheroid`
"""
function ellipsoid_im(ig::MIRT_image_geom, ptype::Symbol ; args...)
	xfov = ig.fovs[1]
	yfov = ig.fovs[2]
	zfov = ig.fovs[3]

	if ptype == :zhu
		params = shepp_logan_3d_parameters(xfov, yfov, zfov, :zhu)
	elseif ptype == :kak
		params = shepp_logan_3d_parameters(xfov, yfov, zfov, :kak)
#=
	elseif ptype == :e3d
		params = shepp_logan_3d_parameters(xfov, yfov, zfov, :e3d)
=#
	elseif ptype == :spheroid
		params = spheroid_params(xfov, yfov, zfov, ig.dx, ig.dy, ig.dz)
	else
		throw("bad phantom symbol $ptype")
	end
	return ellipsoid_im(ig, params ; args...)
end


"""
`phantom = ellipsoid_im(ig ; args...)`

`:zhu` (default) for given image geometry `ig`
"""
function ellipsoid_im(ig::MIRT_image_geom ; args...)
	return ellipsoid_im(ig, :zhu ; args...)
end


"""
`rot3`
"""
function rot3(x, y, z, azim, polar)
	polar != 0 && throw("z (polar) rotation not done")
	xr = cos(azim) * x + sin(azim) * y
	yr = -sin(azim) * x + cos(azim) * y
	zr = z
	return (xr, yr, zr)
end


"""
`shepp_logan_3d_parameters()`

most of these values are unitless 'fractions of field of view'
"""
function shepp_logan_3d_parameters(xfov, yfov, zfov, ptype)
	# parameters from Kak and Slaney, typos ?
	ekak = Float32[
		0		0		0		0.69	0.92	0.9		0	0	2.0
		0		0		0		0.6624	0.874	0.88	0	0	-0.98
		-0.22	0		-0.25	0.41	0.16	0.21	0	0	-0.98
		0.22	0		-0.25	0.31	0.11	0.22	72	0	-0.02
		0		0.1		-0.25	0.046	0.046	0.046	0	0	0.02
		0		0.1		-0.25	0.046	0.046	0.046	0	0	0.02
		-0.8	-0.65	-0.25	0.046	0.023	0.02	0	0	0.01
		0.06	-0.065	 -0.25	0.046	0.023	0.02	90	0	0.01
		0.06	-0.105	 0.625	0.56	0.04	0.1		90	0	0.02
		0		0.1		-0.625	0.056	0.056	0.1		0	0	-0.02
	]

	#parameters from leizhu@standford.edu, says kak&slaney are incorrect
	ezhu = Float32[
	0		0		0		0.69	0.92	0.9		0	0	2.0
	0		-0.0184	0		0.6624	0.874	0.88	0	0	-0.98
	-0.22	0		-0.25	0.41	0.16	0.21	-72	0	-0.02
	0.22	0		-0.25	0.31	0.11	0.22	72	0	-0.02
	0		0.35	-0.25	0.21	0.25	0.35	0	0	0.01
	0		0.1		-0.25	0.046	0.046	0.046	0	0	0.01
	-0.08	-0.605	-0.25	0.046	0.023	0.02	0	0	0.01
	0		-0.1	-0.25	0.046	0.046	0.046	0	0	0.01
	0		-0.605	-0.25	0.023	0.023	0.023	0	0	0.01
	0.06	-0.605	-0.25	0.046	0.023	0.02	-90	0	0.01
	0.06	-0.105	0.0625	0.056	0.04	0.1		-90	0	0.02
	0		0.1		0.625	0.056	0.056	0.1		0	0	-0.02
	]

	if ptype == :zhu
		params = ezhu
	elseif ptype == :kak
		params = ekak
#=
	elseif ptype == :shepp_logan_e3d || ptype == :e3d
		params = e3d[:, [5:7 2:4 8 1]]
=#
	else
		throw("unknown parameter type $ptype")
	end

	params[:,[1,4]] .*= xfov/2
	params[:,[2,5]] .*= yfov/2
	params[:,[3,6]] .*= zfov/2
	return params
end


"""
`spheroid_params()`
"""
function spheroid_params(xfov, yfov, zfov, dx, dy, dz)
#	xfov = nx * dx, number * size
	xradius = (xfov/2) - dx
	yradius = (yfov/2) - dy
	zradius = (zfov/2) - dz
	return [ 0 0 0 xradius yradius zradius 0 0 1 ]
end


"""
`ellipsoid_im()`

show docstring
"""
function ellipsoid_im()
	@doc ellipsoid_im
end


"""
`ellipsoid_im_show()`
"""
function ellipsoid_im_show()
	ig = image_geom(nx=512, nz=64, dz=0.625, fov=500)
	ig = ig.down(8)

	spheroid = ellipsoid_im(ig, :spheroid)
	p1 = jim(spheroid, title="spheroid")

	x = ellipsoid_im(ig, :zhu; hu_scale=1000)
	p2 = jim(x, title="zhu") #clim=[900,1100]

	plot(p1, p2)
end


"""
`ellipsoid_im_test()`
"""
function ellipsoid_im_test()
	ig = image_geom(nx=512, nz=64*2, dz=0.625, fov=500)
	ig = ig.down(8)

	#test different ellipses
	ellipsoid_im(ig)
	ellipsoid_im(ig, :zhu)
	ellipsoid_im(ig, :kak)
	ellipsoid_im(ig, :spheroid, return_params=true)
	ellipsoid_im(ig; checkfov=true)
	ellipsoid_im(ig; showmem=true)
	@test_throws String ellipsoid_im(ig, :bad)
	@test_throws String ellipsoid_im(ig, :spheroid, how=:bad)
	@test_throws String shepp_logan_3d_parameters(0, 0, 0, :bad)
	tmp = [1, 1, 1, 0, 0, 0]
	@test !ellipsoid_im_check_fov(2, 2, 2, [3 0 0 1 1 1 0 0 1], tmp...)
	@test !ellipsoid_im_check_fov(2, 2, 2, [0 3 0 1 1 1 0 0 1], tmp...)
	@test !ellipsoid_im_check_fov(2, 2, 2, [0 0 3 1 1 1 0 0 1], tmp...)

	#ell1 = ellipsoid_im(ig, :zhu; how=:fast) # fast doesn't work
	#ell2 = ellipsoid_im(ig, :zhu; how=:lowmem) # lowmem calls fast - doesn't work
	# todo: test slow vs fast
	true
end


"""
`ellipsoid_im(:test)`

`ellipsoid_im(:show)`

run tests
"""
function ellipsoid_im(test::Symbol)
	if test == :show
		return ellipsoid_im_show()
	end
	test != :test && throw(ArgumentError("test $test"))
	ellipsoid_im()
	ellipsoid_im(:show)
	@test ellipsoid_im_test()
	true
end

#ellipsoid_im(:show)
#ellipsoid_im(:test)
