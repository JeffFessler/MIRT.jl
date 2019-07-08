using Plots

"""
`(sino, pos, ang) = ellipse_sino(sg, ells;
            oversample=1, xscale=1, yscale=1)`

Create sinogram projection of one or more ellipses.
Works for both parallel-beam geometry and for fan-beam geometry.

in
* `sg`                      sinogram geometry object from sino_geom()
* `ells`            [ne 6]  parameters
                            [centx centy radx rady angle_degrees amplitue]
                            (or strum motion object; see ellipse_motion.m)

options
* `oversample::Integer`     oversampling factor for emulating "strips"
                            default 1: just 1 ray per detector element
* `xscale::Integer`         use -1 to flip in x (not recommended)
* `yscale::Integer`         use -1 to flip in y (not recommended)
* `type::Symbol`            `:fan` `:par` `:moj`

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
        type::Symbol=:fan) # what should be the default type

	type = sg.type
        if type == :fan
                (sino, pos, ang) = ellipse_sino_go(ells, [], [],
                                        sg.nb, sg.ds, sg.offset_s, sg.na, sg.orbit, sg.orbit_start,
                                        sg.dso, sg.dod, sg.dfs, sg.sourse_offset, xscale, yscale, oversample, 0)
        elseif type == :par
                (sino, pos, ang) = ellipse_sino_go(ells, [], [],
                                        sg.nb, sg.dr, sg.offset_r, sg.na, sg.orbit, sg.orbit_start,
                                        Inf, 1, 0, sg.source_offset, xscale, yscale, oversample, 0)
        elseif type == :moj
                (sino, pos, ang) = ellipse_sino_go(ells, [], [],
					sg.nb, [], sg.offset_r, sg.na, sg.orbit, sg.orbit_start,
					Inf, 1, 0, sg.source_offset, xscale, yscale, oversample, sg.dx)
	else
		throw("sino geom $type not done")
	end
        return (sino, pos, ang)
end

"""
`ellipse_sino_go()`
"""
function ellipse_sino_go(ells, pos, ang, nb, ds, offset_s, na, orbit, orbit_start,
	dso, dod, dfs, source_offset, xscale, yscale, oversample, mojette)

	if isempty(sng)
		ang = ((orbit_start + [0:na - 1]/na * orbit)) * (pi/180)
	end

	(pos, pos2) = ellipse_sino_pos(pos[:], nb, ds, offset_s, oversample, mojette, ang)

	sino = ellipse_sino_do(ells, pos2, ang[:]', xscale, yscale, dso, dod, dfs, source_offset)

	if oversample != 0
		sino = downsample2(sino, [oversample 1])
	end

	return (sino, pos, ang)
end

"""
`ellipse_sino_pos()`

determine usual and fine "radial" sample positions
"""
function ellipse_sino_pos(pos, nb, ds, offset_s, nover, mojette, ang)

	if isempty(pos)
		if isempty(nb)
			throw("nb required when no pos provided")
		end
		wb = (nb - 1)/2 + offset_s
	else
		if !isempty(nb)
			throw("nb ignored when pos provided")
		end
	end

	if mojette != 0 # tricky mojette radial sampling
		# trick: ray_spacing aka ds comes from dx which is in mojette
		if !isempty(ds)
			throw("ds must be empty for mojette case")
		end
		dt = (abs(mojette) * max(abs(cos(ang))), abs(sin(ang)))

		if !isempty(pos)
			throw("mojette requires empty 'pos'")
		end

		na = max(size(ang))
		pos_coarse = ([0:nb-1] - wb)' * dt[:]' # [nb na]
		if nover > 1
			pos_fine = zeros(nover*nb, na)
			for ia in 1:na
				tmp = [-(nover-1):2:(nover-1)]' / (2 * nover) * dt[ia]
				#pos_fine[:, ia] =
				#matlab: pos_fine(:,ia) = col(outer_sum(tmp, pos_coarse(:,ia)'));

			end
		else
			pos_fine = pos_coarse
		end
	else # ordinary non-mojette sampling
		if isempty(pos)
			pos_coarse = ([0:nb - 1]' - wb) * ds # [nb 1]
		else
			pos_coarse = pos[:] # [nb 1]
			ds = pos[2] - pos[1]
			#=
			if any(abs(diff(pos) / ds - 1) > 1e-10)
				error 'uniform spacing required'
			end =#
		end

		if nover > 1
			# determine fine sampling positions
			# to do: allow different case for trapezoidal rule
			pos_fine = [-(nover - 1):2:(nover - 1)]' / (2*nover) * ds
			#pos_fine = outer_sum(pos_fine, pos_coarse[:]')
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

	nb = size(pos)[1]
	na = max(size(ang))

	# effective radial and angular sample locations in parallel geometry

	if isinf(dso) # type of dso?
		if size(pos)[2] > 1 #mojette
			rads = pos
			angs = repeat(ang, nb, 1)
		else
			(rads, angs) = ndgrid(pos, ang)
		end
	else # fan
		if size(pos)[2] > 1
			throw("mojette fan not supported")
		end
		dis_src_det = dso + dod

		if isinf(dfs) # type of dfs?
			gam = atan(pos / dis_src_det) # gamma
		else # arc detector
			dis_foc_det = dfs + dis_src_det
			ald = pos / dis_foc_det
			gam = atan(dis_foc_det * sin(alf), dis_foc_det * cos(alf) - dfs) # gamma
		end

		rad = dso * sin(gam) + source_offset * cos(gam)
		rads = rad[:, ones(1, na)] # [nb na]
		angs = outer_sum(gam, ang) # [nb na] gamma + beta
	end

	# clear alf gam rad pos ang

	sino = zeros(nb, na)

	cangs = cos.(angs)
	sangs = sin.(angs)

	#loop over ellipses
	#ticker reset
	#if isa(ells, "strum")
		#ne = ells.ne
	#else
	ne = size(ells)[1]
	if size(ells)[2] != 6
		throw("6 parameters per ellipse")
	end
	#end
	for ie in 1:ne
		#ticker(mfilename, ie, ne)
		if isa(ells, "strum")
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
		eang = ell[:, :, 5] .* (pi/180)
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
		sino += scale ./ rp2 .* sqrt(max(rp2 - dis2, 0))
	end
	return sino
end


# ellipse sino old

# ellipse sino orig

"""
`ellipse_sino_test()`

internal test routine: standard sampling
"""
function ellipse_sino_test()
	down = 4
	ig = image_geom(nx=512, ny=504, dx=1)
	ell = [40 70 50 150 20 10]
	(xtrue, ell) = ellipse_im(ig, ell, oversample = 4)

	gf = sino_geom(:fan, nb = 888, na = 984, ds = 1.0, offset_s = 0.75, orbit = 360,
			orbit_start = 0, dsd = 949, dod = 408, down=down)
	# dfs = Inf, source_offset = 0.7, flat fan, not working
	gp = sino_geom(:par, nb = 888, na = 984, dr = 0.5, offset_r = 0.25, orbit = gf.orbit,
			orbit_start = gf.orbit_start, down = down)

	oversample = 8
	sino_mf = ellipse_sino(gf, ell, oversample=oversample) # fan
	sino_mp = ellipse_sino(gf, ell, oversample=1) # parallel

	p1 = jim(sino_mf, title="fan")

	plot(p1)
end

ellipse_sino_test()
