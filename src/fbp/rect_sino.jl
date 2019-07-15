#=
rect_sino.jl
2019-07-12, Helena H
2019-07-13, Jeff Fessler, refactor to use sg.grid
=#

export rect_sino

using Plots: plot
#using MIRT: sino_geom, MIRT_sino_geom, downsample2


"""
`sino = rect_sino(sg, rects ; oversample=1, xscale=1, yscale=1)`

Create sinogram projection of one or more rectangles.
Works for any sinogram geometry.

in
- `sg::MIRT_sino_geom,`		sinogram geometry object from `sino_geom()`
- `rects::Matrix`			`[ne 6]` rectangle parameters
							`[centx centy widthx widthy angle_degrees value]`

options
- `oversample::Integer`	oversampling factor for emulating "strips"
						default 1: just 1 ray per detector element
- `xscale::Integer`		use -1 to flip in x (not recommended); default 1
- `yscale::Integer`		use -1 to flip in y (not recommended); default 1

out
- `sino`		`[nb na]` sinogram

To get the sample locations, use `(pos,ang) = sg.grid`
"""
function rect_sino(sg::MIRT_sino_geom, rects::AbstractMatrix{<:Real} ;
		oversample::Integer=1, kwargs...)

	sg = sg.over(oversample)
	(rg, ϕg) = sg.grid
	sino = rect_sino(rg, ϕg, rects; kwargs...)
	if oversample > 1
		sino = downsample2(sino, [oversample 1])
	end
	return sino
end


"""
`sino = rect_sino(rg::AbstractArray{<:Real}, ϕg::AbstractArray{<:Real},
		rects ; oversample=1, xscale=1, yscale=1)`

Create sinogram projection of one or more rectangles,
for arbitrary radial/angular sampling grid locations `(rg, ϕg)`

in
- `rg::AbstractArray{<:Real}`	radial sampling locations
- `ϕg::AbstractArray{<:Real}`	angular sampling locations (radians)
- `rects::Matrix`			`[n 6]` rectangle parameters
							`[centx centy widthx widthy angle_degrees value]`

options
- `xscale::Integer`		use -1 to flip in x (not recommended); default 1
- `yscale::Integer`		use -1 to flip in y (not recommended); default 1

out
- `sino::AbstractArray{Float32}` same size as `rg` and `ϕg`
"""
function rect_sino(rg::AbstractArray{<:Real}, ϕg::AbstractArray{<:Real},
		rects::AbstractMatrix{<:Real} ;
		xscale::Integer=1,
		yscale::Integer=1,
	)

	size(rects,2) != 6 && throw("6 parameters per rectangle")
	size(rg) != size(ϕg) && throw("rg and ϕg size mismatch")

	sino = zeros(Float32, size(rg))

	cangs = cos.(ϕg)
	sangs = sin.(ϕg)

	#loop over rects
	#ticker reset
	ne = size(rects, 1)
	for ie in 1:ne
		#ticker(mfilename, ie, ne)
		rect = rects[ie, :]

		cx = rect[1] * xscale
		cy = rect[2] * yscale
		wx = rect[3]
		wy = rect[4]
		(wx <= 0) || (wy <= 0) && throw("need positive rectangle sizes")
		eang = deg2rad(rect[5])
		val = rect[6]

		if yscale == -1
			eang = -eang
		end
		if xscale == -1
			eang = pi - eang
		end

		cos_rot = cangs .* cos(eang) + sangs .* sin(eang)
		sin_rot = sangs * cos(eang) - cangs * sin(eang)
		rp = sqrt.((wx * cos_rot).^2 + (wy * sin_rot).^2) # projected radius

		sp = cx .* cangs + cy .* sangs # radial shift
		dis = (rg - sp) ./ rp # scaled distance from center

		# projection angle after affine scaling and rotate
		abs_cos_ang_pi = wx * abs.(cos_rot) ./ rp
		abs_sin_ang_pi = wy * abs.(sin_rot) ./ rp

		# break points of the trapezoid
		len = 1 ./ max.(abs_cos_ang_pi, abs_sin_ang_pi)
		dmax = (abs_cos_ang_pi + abs_sin_ang_pi) / 2
		dbreak = abs.(abs_cos_ang_pi - abs_sin_ang_pi) / 2
		dmax_break = dmax - dbreak
		scale = val * wx * wy ./ rp .* len #lmax

		sino += scale .* trapezoid.(dis, -dmax, -dbreak, dbreak, dmax)
	end

	return sino
end


"""
`trapezoid(t::Real, t1, t2, t3, t4)`
"""
function trapezoid(t::Real, t1::Real, t2::Real, t3::Real, t4::Real)
	if t1 < t < t2
		return (t - t1)/(t2 - t1)
	elseif t2 <= t <= t3
		return 1
	elseif t3 < t < t4
		return (t4 - t)/(t4 - t3)
	else
		return 0
	end
end


"""
`rect_sino()`

show doc strings
"""
function rect_sino()
	@doc rect_sino
end


"""
`rect_sino_test()`

internal test routine: standard sampling
"""
function rect_sino_test()
	rect = [60 10.5 38 27 0 1]
	sg = sino_geom(:ge1, down=8)
	rect_sino(sg, rect; xscale=-1, yscale=-1) # test scale
	t = LinRange(-2, 5, 101)
	f = trapezoid.(t, 1, 2, 3, 4.)
#	plot(t,f)
	true
end


"""
`rect_sino_show()`
show examples
"""
function rect_sino_show( ;
		down::Int = 4,
		rect::Matrix = [70 20 70 30 0 1],
		orbit::Real = 360,
		na::Int = 400,
		oversample::Int = 2,
	)

#=
	ig = image_geom(nx=512, ny=504, fov=500)
	ig = ig.down(down)
	xtrue = rect_im(ig, rect; oversample=4)
=#

	arg = (na=na, down=down, orbit=orbit, offset=0.25)
	geoms = (
		sino_geom(:par, nb = 444, d = 1 ; arg...),
		sino_geom(:fan, nb = 888, d = 1 ; arg..., dsd = 949, dod = 408),
		sino_geom(:fan, nb = 888, d = 1 ; arg...,
			dsd = 949, dod = 408, dfs = Inf), # flat fan
		sino_geom(:moj, nb = 666, d = 1 ; arg...),
	)

	ngeom = length(geoms)
	pl = Array{Plot}(undef, ngeom)

	for ii=1:ngeom
		sg = geoms[ii]
		sino = rect_sino(sg, rect; oversample=oversample)
		dfs = sg.how == :fan ? " dfs=$(sg.dfs)" : ""
		pl[ii] = jim(sg.r, sg.ad, sino, title="$(sg.how)$dfs")
	end
	plot(pl...)
end


"""
`rect_sino(:test)`
self test
"""
function rect_sino(test::Symbol)
	test == :show && return rect_sino_show()
	test != :test && throw(ArgumentError("test $test"))
	rect_sino() # doc
	rect_sino_test()
	rect_sino(:show)
	true
end
