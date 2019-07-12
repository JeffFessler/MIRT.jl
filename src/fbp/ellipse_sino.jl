using Plots
using Test: @test_throws

export ellipse_sino

"""
`(sino, pos, ang) = ellipse_sino(sg, ells;
            oversample=1, xscale=1, yscale=1)`

Create sinogram projection of one or more ellipses.
Works for both parallel-beam geometry and for fan-beam geometry.

in
* `sg`                      sinogram geometry object from sino_geom()
* `ells`            		`[ne 6]`  parameters
                            [centx centy radx rady angle_degrees amplitue]
                            (or strum motion object; see ellipse_motion.m)

options
* `oversample::Integer`     oversampling factor for emulating "strips"
                            default 1: just 1 ray per detector element
* `xscale::Integer`         use -1 to flip in x (not recommended)
* `yscale::Integer`         use -1 to flip in y (not recommended)


out
* `sino`            `[nb na]` sinogram
* `pos`             `[nb 1]`  radial sample positions (or [nb na] if mojette)
* `ang`             `[na]`    angular sample locations
"""
function ellipse_sino(sg::MIRT_sino_geom,
        ells::AbstractArray{<:Real,2};
        oversample::Integer=1,
        xscale::Integer=1,
        yscale::Integer=1,
        )

		how = sg.how
        if how == :fan
                (sino, pos, ang) = ellipse_sino_go(ells, sg.nb, sg.ds, sg.offset, sg.na, sg.orbit, sg.orbit_start,
                                        sg.dso, sg.dod, sg.dfs, sg.source_offset, xscale, yscale, oversample, 0)
        elseif how == :par
                (sino, pos, ang) = ellipse_sino_go(ells, sg.nb, sg.dr, sg.offset, sg.na, sg.orbit, sg.orbit_start,
                                        Inf, 1, 0, sg.source_offset, xscale, yscale, oversample, 0)
        elseif how == :moj
                (sino, pos, ang) = ellipse_sino_go(ells, sg.nb, sg.d, sg.offset, sg.na, sg.orbit, sg.orbit_start,
					Inf, 1, 0, sg.source_offset, xscale, yscale, oversample, 0)
	else
		throw("sino geom $how not done")
	end
        return (sino, pos, ang)
end

"""
`ellipse_sino_go()`
"""
function ellipse_sino_go(ells, nb, ds, offset_s, na, orbit, orbit_start,
	dso, dod, dfs, source_offset, xscale, yscale, oversample, mojette)

	ang = ((orbit_start .+ (0:(na - 1))' / na * orbit)) * (pi/180)
	(pos, pos2) = ellipse_sino_pos(pos[:], nb, ds, offset_s, oversample, mojette, ang)
	sino = ellipse_sino_do(ells, pos2, ang[:]', xscale, yscale, dso, dod, dfs, source_offset)

	if oversample != 1
		sino = downsample2(sino, [oversample 1])
	end

	return (sino, pos, ang)
end

"""
`ellipse_sino_pos()`

determine usual and fine "radial" sample positions
"""
function ellipse_sino_pos(pos, nb, ds, offset_s, nover, mojette, ang)
	wb = (nb - 1)/2 + offset_s
	if mojette != 0 # tricky mojette radial sampling
		# trick: ray_spacing aka ds comes from dx which is in mojette
		!isempty(ds) && throw("ds must be empty for mojette case")
		dt = (abs(mojette) * max(abs(cos(ang))), abs(sin(ang)))

		!isempty(pos) && throw("mojette requires empty 'pos'")

		na = max(size(ang))
		pos_coarse = ((0:(nb-1)) - wb) * dt[:]' # [nb na]
		if nover > 1
			pos_fine = zeros(nover*nb, na)
			for ia in 1:na
				tmp = (-(nover-1):2:(nover-1)) / (2 * nover) * dt[ia]
				pos_fine[:, ia] = (tmp .+ pos_coarse[:,ia]')[:]

			end
		else
			pos_fine = pos_coarse
		end
	else # ordinary non-mojette sampling
		pos_coarse = ((0:(nb - 1)) .- wb) .* ds' # [nb]
		if nover > 1
			# determine fine sampling positions
			# to do: allow different case for trapezoidal rule
			pos_fine = (-(nover - 1):2:(nover - 1)) / (2*nover) * ds
			pos_fine = pos_fine .+ pos_coarse'
			pos_fine = pos_fine[:]
		else
			pos_fine = pos_coarse
		end
	end

	return (pos_coarse, pos_fine)
end


"""
`ellipse_sino_do()`

analytical line-integral projections of ellipse (fan-beam in general)
"""
function ellipse_sino_do(ells, pos, ang, xscale, yscale, dso, dod, dfs, source_offset)
	nb = size(pos, 1)
	na = maximum(size(ang))
	# effective radial and angular sample locations in parallel geometry

	if isinf(dso) # type of dso?
		if size(pos, 2) > 1 #mojette
			rads = pos
			angs = repeat(ang', nb, 1)
		else
			(rads, angs) = ndgrid(pos, ang')
		end

	else # fan
		size(pos, 2) > 1 && throw("mojette fan not supported")

		dis_src_det = dso + dod

		if isinf(dfs) # type of dfs?
			gam = atan(pos / dis_src_det) # gamma
		else # arc detector
			dis_foc_det = dfs + dis_src_det
			alf = pos / dis_foc_det
			gam = atan.(dis_foc_det * sin.(alf), dis_foc_det * cos.(alf) .- dfs) # gamma
		end

		rad = dso * sin.(gam) + source_offset * cos.(gam)
		rads = repeat(rad, 1, na) # [nb, na]
		angs = gam .+ ang  # [nb na] gamma + beta
	end

	# clear alf gam rad pos ang

	sino = zeros(nb, na)

	cangs = cos.(angs)
	sangs = sin.(angs)

	#loop over ellipses
	#ticker reset

	ne = size(ells)[1]
	size(ells)[2] != 6 && throw("6 parameters per ellipse")
	for ie in 1:ne
		#ticker(mfilename, ie, ne)
		if isa(ells, Type)
			ell = ells.ell[ie, na] # [na 6]
			ell = reshape(ell, 1, na, 6)
			ell = repeat(ell, nb, 1, 1)
		else
			ell = ells[ie, :]
			ell = reshape(ell, 1, 1, 6) # is this 3 dimensional
		end

		cx = xscale * ell[:, :, 1]
		rx = ell[:, :, 3]
		cy = yscale * ell[:, :, 2]
		ry = ell[:, :, 4]
		eang = ell[5] * (pi/180)
		val = ell[:, :, 6]

		if yscale == -1
			eang = -eang
		end
		if xscale == -1
			eang = pi - eang
		end
		scale = 2 * val .* rx .* ry

		# square of projected radius:
		rp2 = (rx .* (cangs .* cos(eang) + sangs .* sin(eang))).^2 + (ry .* (sangs .* cos(eang) - cangs .* sin(eang))).^2
		sp = cx .* cangs + cy .* sangs # radial shift
		dis2 = (rads - sp).^2 # square of distances from center
		sino += scale ./ rp2 .* sqrt.(max.(rp2 - dis2, 0))
	end
	return sino
end

"""
`ellipse_sino()`

shows doc strings
"""
function ellipse_sino()
	@doc ellipse_sino
end

"""
`ellipse_sino_test()`

internal test routine: standard sampling
"""
function ellipse_sino_test()
	down = 4
	ig = image_geom(nx=512, ny=504, dx=1)
	ell = [
	40 70 50 150 20 10
	]
	#(xtrue, ell) = ellipse_im(ig, ell; oversample=4)


	gf = sino_geom(:fan, nb = 888, na = 984, d = 1.0, orbit = 360,
			orbit_start = 0, offset = 0.75, dsd = 949, dod = 408, down=down)
	# dfs = Inf, source_offset = 0.7, flat fan, not working
	gp = sino_geom(:par, nb = 888, na = 984, down=down, d = 0.5, orbit = gf.orbit,
			orbit_start = gf.orbit_start, offset = 0.25)
	gm = sino_geom(:moj, nb = 888, na = 984, down=down, d = 0.5, orbit = gf.orbit,
			orbit_start = gf.orbit_start, offset = 0.25)
	sino = sino_geom(:bad)
	@test_throws String ellipse_sino(sino, ell)

	oversample = 8

	sino_mf = ellipse_sino(gf, ell; oversample=oversample) # fan
	sino_mp = ellipse_sino(gp, ell; oversample=1) # parallel
	sino_m = ellipse_sino(gm, ell; oversample=1) # mojette
	e1 = ellipse_sino(gf, ell; oversample=1, xscale=-1, yscale=-1)
	e2 = ellipse_sino(gp, ell; oversample=1, xscale=-1, yscale=-1)
	e3 = ellipse_sino(gm, ell; oversample=1, xscale=-1, yscale=-1)
end

"""
`ellipse_sino_show()`
"""
function ellipse_sino_show()
	down = 4
	ig = image_geom(nx=512, ny=504, dx=1)
	ell = [
	40 70 50 150 20 10
	]

	#(xtrue, ell) = ellipse_im(ig, ell; oversample=4)

	gf = sino_geom(:fan, nb = 888, na = 984, d = 1.0, orbit = 360,
			orbit_start = 0, offset = 0.75, dsd = 949, dod = 408, down=down)
	# dfs = Inf, source_offset = 0.7, flat fan, not working
	gp = sino_geom(:par, nb = 888, na = 984, down=down, d = 0.5, orbit = gf.orbit,
			orbit_start = gf.orbit_start, offset = 0.25)

	oversample = 8

	sino_mf = ellipse_sino(gf, ell; oversample=oversample)[1] # fan
	sino_mp = ellipse_sino(gp, ell; oversample=1)[1] # parallel
	p1 = jim(sino_mf, title="fan")
	p2 = jim(sino_mp, title="parallel")

	plot(p1, p2)
end


"""
`ellipse_sino(:test)`

`ellipse_sino(:show)`
"""
function ellipse_sino(test::Symbol)
	if test == :show
		return ellipse_sino_show()
	end
	ellipse_sino()

	ellipse_sino_test()
	ellipse_sino(:show)
	true
end

#ellipse_sino(:show)
