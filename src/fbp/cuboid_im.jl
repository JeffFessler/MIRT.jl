#=
cuboid_im.jl
2019-07-01 Helena
=#

export cuboid_im

using Plots
using Test: @test, @test_throws
#using MIRT: jim, MIRT_image_geom, downsample3, rotate2d


"""
`phantom = cuboid_im(ig, params;
	oversample=1, hu_scale=1, how=:auto, return_params=false)`

generate cuboid phantom image from parameters:
	`[x_center y_center z_center x_diameter y_diameter z_diameter
		xy_angle_degrees z_angle_degrees amplitude]`

in
- `ig::MIRT_image_geom`		`image_geom()` object
- `params`					`[N 9]`  cuboid parameters.
							note: "diameter" not "radius"

options
- `oversample::Integer` 	 oversampling factor (for partial volume)
- `how::Symbol`				`:sample` use samples
							`:lowmem1` one slice per time
							`:exact` perfect partial volume if angle* = 0
							default: `:exact` if non rotated, else `:sample`
- `return_params::Bool`	if true, return both phantom and params

out
- `phantom`		`[nx ny nz]` image
- `params`		`[N 9]` cuboid parameters (only return if return_params=true)
"""
function cuboid_im(ig::MIRT_image_geom, params::AbstractMatrix{<:Real} ;
		oversample::Integer = 1,
		how::Symbol = :auto,
		return_params::Bool = false)

	size(params,2) != 9 && throw("bad cuboid parameter vector size")
	any(params[:,4:6] .< 0) && throw("cuboid `diameters` must be nonnegative")

	if oversample > 1
		ig = ig.over(oversample)
	end
	args = (ig.nx, ig.ny, ig.nz, params,
		ig.dx, ig.dy, ig.dz, ig.offset_x, ig.offset_y, ig.offset_z)

	phantom = zeros(Float32, ig.nx, ig.ny, ig.nz)

	if how == :auto
		rotated = params[:, 7:8] .!= 0
		if any(rotated)
			how = :sample
		else
			how = :exact
		end
	end

	if how == :exact
		phantom += cuboid_im_exact(args...)
	elseif how == :lowmem1
		phantom += cuboid_im_lowmem1(args...)
	elseif how == :sample
		phantom += cuboid_im_sample(args...)
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
`cuboid_im_sample()``

:sample
"""
function cuboid_im_sample(nx, ny, nz, params,
		dx, dy, dz, offset_x, offset_y, offset_z)
	phantom = zeros(Float32, nx, ny, nz)

	wx = (nx-1)/2 + offset_x
	wy = (ny-1)/2 + offset_y
	wz = (nz-1)/2 + offset_z
	x = ((0:nx-1) .- wx) * dx
	y = ((0:ny-1) .- wy) * dy
	z = ((0:nz-1) .- wz) * dz
	(xx, yy, zz) = ndgrid(x, y ,z)

	# ticker reset
	ne = size(params)[1]
	for ie in 1:ne
		# ticker(ie, ne)

		par = params[ie, :]
		cx = par[1]
		cy = par[2]
		cz = par[3]
		rx = par[4]
		ry = par[5]
		rz = par[6]

		theta = deg2rad(par[7])
		phi = deg2rad(par[8])
		phi != 0 && throw("phi rotation not done")

		(x, y) = rotate2d(xx .- cx, yy .- cy, theta)
		z = zz .- cz

		# rx, ry, rz are 'diameters' not 'radius'
		tmp = (abs.(x / rx) .<= 1/2) .& (abs.(y / ry) .<= 1/2) .& (abs.(z / rz) .<= 1/2)
		value = Float32(par[9])
		phantom += value * tmp
	end

	return phantom
end


"""
`cuboid_im_exact()``

:exact

non-rotated cuboid -- rotation not done
"""
function cuboid_im_exact(nx, ny, nz, params, dx, dy, dz,
		offset_x, offset_y, offset_z)
#	over != 1 && throw("oversample=$over for exact!?")
	phantom = zeros(Float32, nx, ny, nz)

	wx = (nx-1)/2 + offset_x
	wy = (ny-1)/2 + offset_y
	wz = (nz-1)/2 + offset_z
	x = ((0:nx-1) .- wx) * dx
	y = ((0:ny-1) .- wy) * dy
	z = ((0:nz-1) .- wz) * dz
	(xx, yy, zz) = ndgrid(x, y, z)

	# length of intersection of [a, b] with [c, d]
	fun1 = (a, b, c, d) -> max(min(d,b) - max(a,c), 0)
	fun2 = (x, h) -> fun1(x - 0.5, x + 0.5, -h, h)

	adx = abs(dx)
	ady = abs(dy)
	adz = abs(dz)

	#ticker reset
	ne = size(params)[1]
	for ie in 1:ne
		#ticker(ie, ne)

		par = params[ie, :]
		cx = par[1]
		cy = par[2]
		cz = par[3]
		rx = par[4]
		ry = par[5]
		rz = par[6]

		(par[7] != 0 || par[8] != 0) && throw("rotation not done")

		# rx, ry, rz are "diameters" not "radius"
		x = (xx .- cx) / adx
		hx2 = rx / adx / 2 # half width in pixels
		y = (yy .- cy) / ady
		hy2 = ry / ady / 2
		z = (zz .- cz) / adz
		hz2 = rz / adz / 2

		fx = fun2.(x, hx2)
		fy = fun2.(y, hy2)
		fz = fun2.(z, hz2)
		tmp = (fx .* fy .* fz)
		value = Float32(par[9])
		phantom += value * tmp
	end

	return phantom
end


"""
`cuboid_im_lowmem1(...)`

This version does `:sample` 1 slice at a time to reduce memory
"""
function cuboid_im_lowmem1(nx, ny, nz, params,
		dx, dy, dz, offset_x, offset_y, offset_z)
	phantom = zeros(Float32, nx, ny, nz)
	for iz in 1:nz
		offset_z_new = (nz - 1)/2 + offset_z - (iz-1)
		phantom[:,:,iz] = cuboid_im_sample(nx, ny, 1, params,
			dx, dy, dz, offset_x, offset_y, offset_z_new)
	end
	return phantom
end


"""
`default_cuboid_parameters()`
"""
function default_cuboid_parameters(xfov, yfov, zfov)
	return [
	0	0		0	60	60	2	0	0	1
	20	20		3	10	20	2	0	0	1
	1.4	-0.5	1	2	2.7	3.2	0	0	1
	]
end


"""
`rotated_cuboid_parameters()`
"""
function rotated_cuboid_parameters(xfov, yfov, zfov)
	return [
	0	0	0	60	60	2	45	0	1
	0	0	0	50	50	3	45	0	1
	]
end


#=
"""
`phantom = cuboid_im(nx, dx, params; args...)`

specifying voxel size `dx` and cuboid `params`
"""
function cuboid_im(nx::Integer, dx::Real, params; args...)
	ig = image_geom(nx=nx, dx=1)
	return cuboid_im(ig, params; args...)
end


"""
`phantom = cuboid_im(nx::Integer, params; args...)`

voxel size `dx=1` and cuboid `params`
"""
function cuboid_im(nx::Integer, params; args...)
	return cuboid_im(nx, 1., params; args...)
end


"""
`phantom = cuboid_im(nx::Integer; ny::Integer=nx,dx::Real=1, args...)`

image of size 'nx' by 'ny' (default 'nx') with specified 'dx' (default 1),
defaults to `:default_cuboid`
"""
function cuboid_im(nx::Integer; ny::Integer=nx, dx::Real=1, args...)
	ig = image_geom(nx=nx, ny=ny, dx=dx)
	return cuboid_im(ig, :default; args...)
end


"""
`phantom = cuboid_im(nx::Integer, ny::integer; args...)`

`:default` of size `nx` by `ny`
"""
function cuboid_im(nx::Integer, ny::Integer; args...)
	return cuboid_im(nx, ny=ny, dx=1.; args...)
end
=#


"""
`phantom = cuboid_im(ig, code, args...)`

`code = :default | :rotate` # add more options
"""
function cuboid_im(ig::MIRT_image_geom, params::Symbol; args...)
	fov = ig.fovs
	if params == :default
		params = default_cuboid_parameters(fov, fov, fov)
	elseif params == :rotate
		params = rotated_cuboid_parameters(fov, fov, fov)
	else
		throw("bad phantom symbol $params")
	end
	return cuboid_im(ig, params; args...)
end


"""
`phantom = cuboid_im(ig; args...)`

`:default` default parameters
"""
function cuboid_im(ig::MIRT_image_geom; args...)
	return cuboid_im(ig, :default; args...)
end



"""
`cuboid_im()`

show docstring
"""
function cuboid_im()
	@doc cuboid_im
end


"""
`cuboid_im_show()`
"""
function cuboid_im_show()
	ig = image_geom(nx=2^6, ny=2^6, nz=2^3, dz=1, fov=240)

	x1 = cuboid_im(ig, :default, how=:exact)
	p1 = jim(x1, title="exact")

	x2 = cuboid_im(ig, :default, how=:sample)
	p2 = jim(x2, title="sample")

	x3 = cuboid_im(ig, :default, how=:lowmem1)
	p3 = jim(x3, title="lowmem")

	x4 = cuboid_im(ig, :rotate)
	p4 = jim(x4, title="rotate") # runs sample

	plot(p1, p2, p3, p4)
end



"""
`cuboid_im_test()`
"""
function cuboid_im_test()
	ig = image_geom(nx=2^3, ny=2^3, nz=2^3, dz=1, fov=240)
	fov = 100
	cuboid_im_show()

	@test_throws String cuboid_im(ig, :bad)
	@test_throws String cuboid_im(ig, :rotate, how=:exact)
	@test_throws String cuboid_im(ig, :rotate, how=:bad)
	@test_throws String cuboid_im(ig, [0 0 0 0 0 0 0 1 1])
	@test cuboid_im(ig, :rotate, return_params=true) isa Tuple

	diam = abs.([2*ig.dx 2.7*ig.dy 3.2*ig.dz])
	params = [1.4 -0.5 1 diam 0 0 1]
	phantom_exact = cuboid_im(ig, params, how=:exact)
	dxyz = abs.(ig.dx * ig.dy * ig.dz)
	vol_true = prod(diam)
	vol_phantom = sum(phantom_exact[:]) * dxyz
	isapprox(vol_true, vol_phantom)

	phantom_sample = cuboid_im(ig, :default, how=:sample, oversample=3)
	vol_phantom2 = sum(phantom_sample[:]) * dxyz
	isapprox(vol_phantom2, vol_true)

	xrl = cuboid_im(ig, :rotate, how=:lowmem1)
	xrs = cuboid_im(ig, :rotate, how=:sample)
	@test xrl == xrs

	# test
#	x4 = cuboid_im(ig, :rotate, how=:exact) exact does not support rotation
	x3 = cuboid_im(ig, :default, how=:lowmem1)
	x5 = cuboid_im(ig)
	true
end


"""
`cuboid_im(:test)`

`cuboid_im(:show)`

run tests
"""
function cuboid_im(test::Symbol)
	if test == :show
		return cuboid_im_show()
	end
	test != :test && throw(ArgumentError("test $test"))
	cuboid_im()
	cuboid_im(:show)
	cuboid_im_test()
	true
end

#cuboid_im(:show)
#cuboid_im(:test)
