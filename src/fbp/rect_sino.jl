using Plots

export rect_sino

"""
`(sino, pos, ang) = rect_sino(sg, rects;
            oversample=1, xscale=1, yscale=1)`

Create sinogram projection of one or more rectangles.
Works for both parallel-beam geometry and for fan-beam geometry.

in
* `sg`                  sinogram geometry object from sino_geom()
* `rects`               [ne 6]
                        [centx centy widthx width y angle_degrees value]

options
* `oversample::Integer` oversampling factor for emulating "strips"
                        default 1: just 1 ray per detector element
* `xscale::Integer`     use -1 to flip in x (not recommended)
* `yscale::Integer`     use -1 to flip in y (not recommended)
* `how::Symbol`        `:fan` `:par` `:moj`

out
* `sino`                [nb na] sinogram
* `pos`                 [nb 1]  radial sample positions (or [nb na] is mojette)
* `ang`                 [na]    angular sample locations (in radians)
"""
function rect_sino(sg::MIRT_sino_geom,
    	rects::AbstractArray{<:Real,2};
    	oversample::Integer=1,
    	xscale::Integer=1,
    	yscale::Integer=1)

	size(rects, 2) != 6 && throw("6 parameters per rect")

    if sg.how == :fan
        (sino, pos) = rect_sino_go(rects, sg.nb, sg.ds, sg.offset, sg.na,
                        sg.ar, sg.dso, sg.dod, sg.dfs, sg.source_offset,
                        xscale, yscale, oversample, 0)
    elseif sg.how == :par
        (sino, pos) = rect_sino_go(rects, sg.nb, sg.dr, sg.offset, sg.na,
                        sg.ar, Inf, 1, 0, sg.source_offset, xscale, yscale, oversample, 0)
    elseif sg.how == :moj
        (sino, pos) = rect_sino_go(rects, sg.nb, [], sg.offset, sg.na,
                        sg.ar, Inf, 1, 0, sg.source_offset, xscale, yscale, oversample, sg.d)
    else
        throw("sino geom $(sg.how) not done") # never should happen!
    end

    return (sino, pos)
end


"""
`rect_sino_go()`
"""
function rect_sino_go(rects, nb, ds, offset_s, na, ang, dso, dod, dfs,
    source_offset, xscale, yscale, oversample, mojette)

    (pos, pos2) = rect_sino_pos([], nb, ds, offset_s, oversample, mojette, ang)

    sino = rect_sino_do(rects, pos2, ang[:]', xscale, yscale, dso, dod, dfs, source_offset)

    if oversample != 1
        sino = downsample2(sino, [oversample 1])
    end
    return (sino, pos)
end

"""
`rect_sino_pos()`

determine usual and fine "radial" sample positions
"""
function rect_sino_pos(pos, nb, ds, offset_s, nover, mojette, ang)
    wb = (nb - 1)/2 + offset_s
    if mojette != 0 # tricky mojette radial sampling
        # trick: ray_spacing aka ds comes from dx which is in mojette
        !isempty(ds) && throw("ds must be empty for mojette case")

        dt = abs(mojette) * max.(abs.(cos.(ang)), abs.(sin.(ang)))' # [1, na]

        !isempty(pos) && throw("mojette requires empty 'pos'")

        na = length(ang)
        pos_coarse = ((0:nb - 1) .- wb) * dt[:]' # [nb na]

        if nover > 1
            pos_fine = zeros(nover*nb, na)
            for ia in 1:na
                tmp = (-(over - 1):2:(nover - 1)) / (2*nover) * dt[ia]
                pos_fine[:, ia] = (tmp .+ pos_coarse[:,ia]')[:]
            end
        else
            pos_fine = pos_coarse
        end
    else # ordinary mojette sampling
        pos_coarse = ((0:(nb - 1)) .- wb) .* ds' # [nb]
        if nover > 1
            # determine fine sampling positions
            pos_fine = (-(nover - 1):2:(nover-1)) / (2*nover) * ds
            pos_fine = pos_fine .+ pos_coarse' # [nover nb]
            pos_fine = pos_fine[:]
        else
            pos_fine = pos_coarse
        end
    end
    return (pos_coarse, pos_fine)
end

"""
`rect_sino_do()`

analytical line-integral projections rectangle
"""
function rect_sino_do(rects, pos, ang, xscale, yscale, dso, dod, dfs, source_offset)
	nb = size(pos, 1)
	na = maximum(size(ang))

	# effective radial and angular sample locations in parallel geometry
	if isinf(dso)
		if size(pos, 2) > 1 # mojette
			rads = pos
			angs = repeat(ang, nb, 1) # [nb na]
		else

			(rads, angs) = ndgrid(pos, ang') # [nb na]
		end
	else # fan
		size(pos, 2) > 1 && throw("mojette fan not supported")
		dis_src_det = dso + dod

		if isinf(dfs) # flat detector
			gam = atan(pos / dis_src_det) # gamma
		else
			dis_foc_det = dfs + dis_src_det
			alf = pos / dis_foc_det
			gam = atan.(dis_foc_det * sin.(alf), dis_foc_det * cos.(alf) .- dfs) #gamma
		end

		rad = dso * sin.(gam) + source_offset * cos.(gam)
		rads = repeat(rad, 1, na) # [nb, na]

		angs = gam .+ ang  # [nb na] gamma + beta
	end


	sino = zeros(nb, na)

	cangs = cos.(angs)
	sangs = sin.(angs)

	#loop over rects
	#ticker reset
	ne = size(rects, 1)
	for ie in 1:ne
		#ticker(mfilename, ie, ne)
		rect = rects[ie, :] # [1, 6]

		cx = xscale * rect[1]
		cy = yscale * rect[2]
		wx = rect[3]
		wy = rect[4]
		(wx <= 0) || (wy <= 0) && throw("need positive rectangle sizes")
		eang = rect[5] * (pi/180)
		if yscale == -1
			eang = -eang
		end
		if xscale == -1
			eang = pi - eang
		end
		val = rect[6]

		cos_rot = cangs .* cos(eang) + sangs .* sin(eang)
		sin_rot = sangs * cos(eang) - cangs * sin(eang)
		rp = sqrt.((wx * cos_rot).^2 + (wy * sin_rot).^2) # projected radius

		sp = cx .* cangs + cy .* sangs # radial shift
		dis = (rads - sp) ./ rp # scaled distance from center

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
'div0'
"""
function div0(x, y)
	if y == 0
		return 0
	else
		return x/y
	end
end

"""
"""
function trapezoid(t::Real, t1, t2, t3, t4)
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

shows doc strings
"""
function rect_sino()
	@doc rect_sino
end


"""
`rect_sino_test()`

internal test routine: standard sampling
"""
function rect_sino_test()
	down = 1
	ig = image_geom(nx=512, ny=504, fov=500)
	ig = ig.down(down)
	rect = [[60 10.5 38 27] .* ig.dx 0 1]
	# (xtrue, rect) = rect_im(ig, rect; oversample=3)

	sgf = sino_geom(:fan, nb = 888, na = 984, d = 1.0, orbit = 360,
			orbit_start = 0, offset = 0.75, dsd = 949, dod = 408, down=down) # fan
	sgp = sino_geom(:par, nb = 888, na = 984, down=down, d = 0.5, orbit = sgf.orbit,
			orbit_start = sgf.orbit_start, offset = 0.25) # parallel
	sgm = sino_geom(:moj, nb = 888, na = 984, down=down, d = 0.5, orbit = sgf.orbit,
			orbit_start = sgf.orbit_start, offset = 0.25) # mojette

	sino_f = rect_sino(sgf, rect; oversample=1)
	sino_p = rect_sino(sgp, rect; oversample=1)
	sino_m = rect_sino(sgm, rect; oversample=2)
	r1 = rect_sino(sgf, rect; oversample=1, yscale=-1, xscale=-1)
	r2 = rect_sino(sgp, rect; oversample=1, yscale=-1, xscale=-1)
	r3 = rect_sino(sgm, rect; oversample=1, yscale=-1, xscale=-1)
	t = LinRange(-2, 5, 101)
	f = trapezoid.(t, 1, 2, 3, 4)
	div0(4, 2)
	div0(1, 0)
end

"""
`rect_sino_show`
"""
function rect_sino_show()
	down = 1
	ig = image_geom(nx=512, ny=504, fov=500)
	ig = ig.down(down)
	rect = [[60 10.5 38 27] .* ig.dx 0 1]

	sgf = sino_geom(:fan, nb = 888, na = 984, d = 1.0, orbit = 360,
			orbit_start = 0, offset = 0.75, dsd = 949, dod = 408, down=down) # fan
	sgp = sino_geom(:par, nb = 888, na = 984, down=down, d = 0.5, orbit = sgf.orbit,
			orbit_start = sgf.orbit_start, offset = 0.25) # parallel

	#sg = sino_geom(:ge1, nb=888, na=984, down=down)
	#sino = rect_sino(sg, rect; oversample=10)

	xtrue = rect_im(ig, rect; oversample=3)
	sino_f = rect_sino(sgf, rect; oversample=1)[1]
	sino_p = rect_sino(sgp, rect; oversample=1)[1]

	x1 = jim(sino_f, title="fan")
	x2 = jim(sino_p, title="parallel")
	x3 = jim(xtrue, title="xtrue")
	#x4 = jim(sino, title="ya")

	plot(x1, x2, x3)
end


"""
`rect_sino(:show)`

`rect_sino(:test)`
"""
function rect_sino(test::Symbol)
	if test == :show
		return rect_sino_show()
	end
	rect_sino()
	rect_sino_test()
	rect_sino(:show)
	true
end

#rect_sino(:test)
#rect_sino(:show)
