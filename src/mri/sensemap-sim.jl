# sensemap-sim.jl
# Simulate coil sensitivity maps
# 2019-06-18, Jeff Fessler, University of Michigan

#using MIRT: jim, image_geom

using Elliptic: ellipke
using Test: @test
import Plots # Plot
using Plots: Plot, plot, plot!, gui, contour!, scatter!, quiver, quiver!,
		clibrary, Arrow
	clibrary(:misc)


"""
`function (smap,info) = ir_mri_sensemap_sim(...)`

Simulate 2D or 3D sensitivity maps for sensitivity-encoded MRI based on
grivich:00:tmf http://doi.org/10.1119/1.19461

This code makes maps for multiple coils,
but does not model coupling between coils,
so most likely it is an approximation at best.

option
* `dims::Dims`		image size; default (64, 64)
* `dx::Real`		pixel/voxel dimension; default: 3
* `dy::Real`		pixel/voxel dimension; default: `dx`
* `dz::Real`		""
* `ncoil::Int`		# of coils total; default 4
* `nring::Int`		# of rings of coils; default 1
* `rcoil::Real`		coil radius; default `dx * nx / 2 * 0.50`
* `dz_coil`			ring spacing in z; default `nz*dz/nring`
			(3D geometry is a cylinder)
* `coil_distance::Real`		distance of coil center from isocenter
	for central ring of coils as a multiple of FOVx,
	where `FOVx=nx*dx`; default 1.2
* `orbit::Real`			default 360 [degrees]
* `orbit_start::Union{Real,AbstractVector{<:Real}} = 0` scalar or `nring` [degrees]
* `scale::Symbol`
  + `:none` (default)
  +	`ssos_center` make SSoS of center = 1

out
* `smap	[dims ncoil]`	simulated sensitivity maps (complex!)
* `info::NamedTuple`	geometry information for plots

All length parameters must have same units (e.g., mm or cm)

Matlab notes:
* 2005-6-20, Jeff Fessler and Amanda Funai, University of Michigan
* 2014-08-19 JF more testing, verifying phase is correct
* 2014-09-09 modified for 3D by Mai Le
* 2016-05-03 JF fixes 
"""
function ir_mri_sensemap_sim(;
		dims::Dims = (64,64), # 2D default
		dx::Real = 3,
		dy::Real = dx,
		dz::Real = 3,
		ncoil::Int = 4,
		nring::Int = 1,
		rcoil::Real = dx * dims[1] / 2 * 0.5,
		dz_coil::Real = ((length(dims) == 3) ? dims[3] : 1) * dz / nring,
		coil_distance::Real = 1.2, # multiplies fov/2
		orbit::Real = 360,
		orbit_start::Union{Real,AbstractVector{<:Real}} = 0,
		scale::Symbol = :none, # or :ssos_center
		chat::Bool = false,
	)

	(length(dims) < 2 || length(dims) > 3) && throw("2D or 3D only")
	nx = dims[1]
	ny = dims[2]
	nz = (length(dims) == 3) ? dims[3] : 0

	coils_per_ring = round(Int, ncoil / nring)
	ncoil != nring * coils_per_ring && throw("nring must be divisor of ncoil")

	(smap, info) = ir_mri_sensemap_sim_do(
			nx, ny, nz,
			dx, dy, dz,
			ncoil, coils_per_ring, rcoil, dz_coil,
			orbit, orbit_start, coil_distance, chat)

	scale_center = (nz == 0) ?
		1 / sqrt(sum(abs.(smap[round(Int,nx/2),round(Int,ny/2),:].^2))) :
		1 / sqrt(sum(abs.(smap[round(Int,nx/2),round(Int,ny/2),round(Int,nz/2),:].^2)))

	if scale == :ssos_center
		smap *= Float32(scale_center)
	elseif scale != :none
		throw("scale $scale")
	end

	return (smap, info)
end


"""
`(smap, info) = ir_mri_sensemap_sim_do()`
"""
function ir_mri_sensemap_sim_do(nx, ny, nz,
		dx, dy, dz, ncoil, ncoilpr, rcoil, dz_coil,
		orbit, orbit_start, coil_distance, chat)

	nring = Int(ncoil / ncoilpr)
	rlist = Float32(rcoil) * ones(Float32, ncoilpr, nring) # coil radii

	zerof = (arg...) -> zeros(Float32, arg...)
	plist = zerof(ncoilpr, nring, 3) # position of coil center [x y z]
	nlist = zerof(ncoilpr, nring, 3) # normal vector (inward) from coil center
	olist = zerof(ncoilpr, nring, 3) # unit vector orthogonal to normal vector in x-y
	ulist = zerof(ncoilpr, nring, 3) # upward vector

	if length(orbit_start) == 1
		orbit_start = repeat([orbit_start], nring)
	elseif length(orbit_start) != nring
		throw("bad orbit_start length")
	end

	# cylindrical coil configuration, like abdominal coils
	alist = deg2rad.(orbit) * (0:(ncoilpr-1)) / ncoilpr # coil angles [radians]
	z_ring = ((1:nring) .- (nring+1)/2) * dz_coil
	for ir = 1:nring
		for ic = 1:ncoilpr
			phi = alist[ic] + deg2rad(orbit_start[ir])
			Rad = max(nx/2 * dx, ny/2 * dy) * coil_distance
			plist[ic,ir,:] = [Rad * [cos(phi), sin(phi)]; z_ring[ir]]
			nlist[ic,ir,:] = -[cos(phi), sin(phi), 0*z_ring[ir]] # cylinder
			olist[ic,ir,:] = [-sin(phi), cos(phi), 0]
			ulist[ic,ir,:] .= [0, 0, 1]
		end
	end

	# object coordinates
	x = Float32.(((1:nx) .- (nx+1)/2) * dx)
	y = Float32.(((1:ny) .- (ny+1)/2) * dy)
	z = (nz > 0) ? Float32.(((1:nz) .- (nz+1)/2) * dz) : [0]
	(xx, yy, zz) = ndgrid(x,y,z)

	smap = zeros(ComplexF32, nx, ny, max(nz,1), ncoilpr, nring)
	for ir = 1:nring
		for ic=1:ncoilpr
			# rotate coordinates to correspond to coil orientation
			zr =	(xx .- plist[ic,ir,1]) .* nlist[ic,ir,1] +
					(yy .- plist[ic,ir,2]) .* nlist[ic,ir,2] +
					(zz .- plist[ic,ir,3]) .* nlist[ic,ir,3]
			xr =	xx .* nlist[ic,ir,2] - yy .* nlist[ic,ir,1]
			yr = zz .- plist[ic,ir,3] # translate along object z axis

#=
			# show coordinates
			jim(x, y, xr) #, xlabel x, ylabel y
			jim(x, y, zr)
=#

			# compute sensitivity vectors in coil coordinates
			tmp = ir_mri_smap1.(xr, yr, zr, rlist[ic,ir])
			sx = [p[1] for p in tmp]
			sy = [p[2] for p in tmp]
			sz = [p[3] for p in tmp]

			# coil response depends on tranverse magnetization only?
			# todo: unsure if this should depend on sy and ulist in 3D
				bx = sz * nlist[ic,ir,1] + sx * olist[ic,ir,1]
				by = sz * nlist[ic,ir,2] + sx * olist[ic,ir,2]
			#	bz = sz * nlist[ic,ir,3] + sx * olist[ic,ir,3]
			smap[:,:,:,ic,ir] = complex.(bx, by)

			if chat && nz == 0 # see field components
				tmp = sqrt.(sx.^2 + sz.^2)
				(xx,yy) = ndgrid(x,y)
				plot(jim(x, y, sx, "sx"),
					jim(x, y, sy, "sy"),
					jim(x, y, sz, "sz"),
					quiver(xx[:], yy[:], title="(sx,sy)",
						quiver=((sx./tmp)[:], (sz./tmp)[:])),
					)
				prompt()
			end

			if chat && nz == 0 # see final field components vs phase
				(xx,yy) = ndgrid(x,y)
				bb = sqrt.(bx.^2 + by.^2)
				jim(x, y, angle.(smap[:,:,1,ic,1]), "phase", color=:hsv)
				quiver!(xx[:], yy[:], title="(bx,by)",
						quiver=((bx./bb)[:], (by./bb)[:]))
				prompt()
			end
		end
	end

	smap *= rlist[1] / Float32(2*pi) # trick: scale so maximum is near unity
	smap = nz == 0 ?  reshape(smap, nx, ny, ncoil) :
		reshape(smap, nx, ny, nz, ncoil)

	info = (x=x, y=y, z=z, dx=dx, dy=dy, dz=dz,
		nlist=nlist, plist=plist, rlist=rlist, olist=olist, ulist=ulist,
		nring=nring, ncoilpr=ncoilpr, rcoil=rcoil)
	return smap, info
end


"""
`ir_mri_sensemap_sim_show3()`
shows coil geometry but not the 3D smap
"""
function ir_mri_sensemap_sim_show3(smap, x, y, z, dx, dy, dz,
	nlist, plist, rlist, olist, ulist, nring, ncoilpr, rcoil)

	pcolor = ['c', 'g', 'r']
	pcolor = i -> pcolor[1+rem(i,3)]

	ir_plot3_cube(x,y,z)
	scatter!(plist[:,:,1], plist[:,:,2], plist[:,:,3], label="") # coil centers

#=
	# coil normals - todo: 3d quiver not yet working
	tmp1 = reshape(plist, :, 3)
	tmp2 = reshape(nlist, :, 3)
	quiver!(tmp1[:,1], tmp1[:,2], tmp1[:,3],
		quiver=(tmp2[:,1], tmp2[:,2], tmp2[:,3]), label="3d")
	prompt()
=#

	if true # coils
		for ir = 1:nring
			for ic = 1:ncoilpr
				tmp = LinRange(0, 2*pi, 50)
				tmp = cos.(tmp) * olist[ic,ir,:]' + sin.(tmp) * ulist[ic,ir,:]'
				tmp = repeat(plist[ic,ir,:]', size(tmp,1), 1) + rcoil * tmp
				plot!(tmp[:,1], tmp[:,2], tmp[:,3], label="")
			#	patch(tmp[:,1], tmp[:,2], tmp[:,3], pcolor(ir),
			#		edgecolor=:none, facealpha=0.5)
			end
		end
	end
	plot!()
end


"""
`ir_plot3_cube(x,y,z)`
"""
function ir_plot3_cube(x,y,z)
	x1 = x[1]
	x2 = x[end]
	y1 = y[1]
	y2 = y[end]
	z1 = z[1]
	z2 = z[end]
	x = [x1,x2,x2,x1,x1,x1,x2,x2,x1,x1]
	y = [y1,y1,y2,y2,y1,y1,y1,y2,y2,y1]
	z = [z1,z1,z1,z1,z1,z2,z2,z2,z2,z2]
	plot(x,y,z, label="", xlabel="x", ylabel="y", zlabel="z")
end


"""
`ir_mri_sensemap_sim_show2()`
returns plot
"""
function ir_mri_sensemap_sim_show2(smap, x, y, dx, dy, nlist, plist, rlist)
	if ndims(smap) == 3
		(nx,ny,ncoil) = size(smap)
	elseif ndims(smap) == 4
		(nx,ny,nz,ncoil) = size(smap)
	else
		throw("unknown ndims(smap) = $(ndims(smap))")
	end

	pl = Array{Plot}(undef, ncoil, 3)
	clim = [0, maximum(abs.(smap))]
	xmax = maximum(abs.([x[:]; y[:]; plist[:,:,[1,2]][:]]))
	ymax = xmax
	for ic=1:ncoil
		tmp = smap[:,:,ic]
		p = jim(x, y, abs.(tmp), clim=clim, "Magnitude $ic")
		plot!(p, xlim=[-1,1]*1.1*xmax, xtick=(-1:1) * nx/2 * dx)
		plot!(p, ylim=[-1,1]*1.1*xmax, ytick=(-1:1) * nx/2 * dy)

	 	scatter!(p, [0], [0], marker=:o, label="", color=:green) # center
	 	scatter!(p, plist[:,:,1][:], plist[:,:,2][:],
			marker=:o, label="", color=:blue) # coil location
		xdir = nlist[ic,1,2]
		ydir = nlist[ic,1,1]
		r = rlist[ic,1]
		plot!(p, plist[ic,1,1].+r*xdir*[-1,1], plist[ic,1,2].+r*ydir*[1,-1],
			label="", color=:blue, linewidth=3) # coil
		pl[ic,1] = p

		ph = angle.(tmp) # show raw phase (understandable with hsv colormap)
		p = jim(x, y, ph, "Phase", clim=[-pi,pi], color=:hsv)
		plot!(p, xlim=[-1,1]*1.1*xmax, xtick=(-1:1) * nx/2 * dx)
		plot!(p, ylim=[-1,1]*1.1*xmax, ytick=(-1:1) * nx/2 * dy)
		pl[ic,2] = p
	end

	ssos = sqrt.(sum(abs.(smap).^2, dims=ndims(smap)))
	ssos = ssos / ssos[Int(end/2),Int(end/2)]

	p = jim(x, y, ssos, "SSoS (normalized)",
		#	xlim=[-1,1]*1.1*xmax,
			xtick=(-1:1) * nx/2 * dx,
			ytick=(-1:1) * ny/2 * dy)
	pl[1,3] = p

	if true # quiver plot for 1st coil
		bx = real(smap[:,:,1])
		by = imag(smap[:,:,1])
		(xx,yy) = ndgrid(x,y)
		p = quiver(xx[:], yy[:], quiver=(bx[:], by[:]),
				aspect_ratio = 1,
				arrow = Arrow(:simple,:head, 0.01, 0.01), # no effect!?
				title = "Field pattern in x-y plane")
		pl[2,3] = p
	end

	for ic=3:ncoil
		pl[ic,3] = plot(xaxis=:off,yaxis=:off,grid=:off) # kludge
	end

	return plot(pl[:]...)
end


"""
`ir_mri_smap_r(r, z)`
function for testing near 0
"""
function ir_mri_smap_r(r, z)
	M = 4 * r / ((1 + r)^2 + z^2) # = k^2, see ellipke
	(K,E) = ellipke(M)
	return 2 * z / r * ((1 + r)^2 + z^2)^(-0.5) *
		((1 + r^2 + z^2) / ((1 - r)^2 + z^2) * E - K)
end


"""
`ir_mri_smap1()`
based on grivich:00:tmf

for a circular coil in "x-y plane" of radius "a"

Note that coil x-y plane is not same as object x-y plane!

Returns `(i,j,k)` components of B vector for each `(x,y,z)` location
"""
function ir_mri_smap1(x, y, z, a)
	x = x / a # normalized units
	y = y / a
	z = z / a
	r = sqrt(x^2 + y^2)
	M = 4 * r / ((1 + r)^2 + z^2) # = k^2, see ellipke
	(K,E) = ellipke(M)

	# the following is B_z in eqn (18) in grivich:00:tmf
	# and same as eqn [10] in wang:00:dop to within constant scale factor
	smap_z = 2 * ((1 + r)^2 + z^2)^(-0.5) *
		(K + (1 - r^2 - z^2) / ((1 - r)^2 + z^2) * E)
	smap_z /= a

#=
	if any(r == 0) # test code to explore when r is near 0
		r0 = LinRange(0,5e-7,101)
		z0 = 0.4
		t0 = ir_mri_smap_r.(r0, z0)
		slope = 3*pi * z0 / ((1+z0^2)^2.5)
		plot(r0, t0); plot!(r0, slope * r0)
	end
=#

	# the following is B_r in eqn (17) in grivich:00:tmf
	smap_r = 2 * z / r * ((1+r)^2 + z^2)^(-0.5) *
		((1 + r^2 + z^2) / ((1-r)^2 + z^2) * E - K)

	if abs(r) < 1e-6
		smap_r = 3 * pi * z / ((1 + z^2).^2.5) * r
	end
	smap_r /= a

	(isnan(smap_r) || isnan(smap_z)) && throw("nan")

	smap_x = smap_r * (r == 0 ? 0 : x / r)
	smap_y = smap_r * (r == 0 ? 0 : y / r)

#	phi = atan2(y, x)
#	smap_x = smap_r * cos(phi)
#	smap_y = smap_r * sin(phi)

	return Float32(smap_x), Float32(smap_y), Float32(smap_z)
end


"""
`ellipke(x::AbstractArray{<:Real})`
"""
function ellipke_(x::AbstractArray{<:Real})
	k = similar(x)
	e = similar(x)
	for i=1:length(x)
		(k[i], e[i]) = ellipke(x[i])
	end
	return (k,e)
end


"""
`ir_mri_sensemap_sim_test0()`
show ellipke
"""
function ir_mri_sensemap_sim_test0()
	ir_mri_smap_r.(5e-7, 0.4)
	m = LinRange(0,1,101)
	(k,e) = ellipke_(m)
	plot(m, k, label="k"); plot!(m, e, label="e")
	plot!(yaxis=[0,3*pi/2])
	plot!(ytick=((0:3)*pi/2, ["0", "\\pi/2", "\\pi", "3\\pi/2"]))
end


"""
`ir_mri_sensemap_sim_test1()`
test ir_mri_smap1 routine, cf Fig. 4 of grivich:00:tmf
"""
function ir_mri_sensemap_sim_test1()
	a = 1
	x = LinRange{Float32}(-2,2,99)
	y = LinRange{Float32}(-2,2,97)
#	zlist = Float32[0.001, 0.1, 0.2, 0.5, 1.0]
	zlist = Float32[0.1, 0.2, 0.5, 1.0]
#	(xx,yy,zz) = ndgrid(x, y, zlist)
	tmp = Iterators.product(x, y, zlist)
	tmp = [ir_mri_smap1(i[1], i[2], i[3], a) for i in tmp]
	smap_x = [p[1] for p in tmp]
	smap_y = [p[2] for p in tmp]
	smap_z = [p[3] for p in tmp]
#	(smap_x, smap_y, smap_z) = ir_mri_smap1.(xx, yy, zz, a)
	smap_b = @. sqrt(smap_x^2 + smap_y^2)
#	return plot(ir_mri_sensemap_sim_test1_show(smap_x, x, y, zlist, "x")...)
	return plot(layout=(4,length(zlist)),
			ir_mri_sensemap_sim_test1_show(smap_x, x, y, zlist, "x")...,
			ir_mri_sensemap_sim_test1_show(smap_y, x, y, zlist, "y")...,
			ir_mri_sensemap_sim_test1_show(smap_z, x, y, zlist, "z")...,
			ir_mri_sensemap_sim_test1_show(smap_b, x, y, zlist, "b")...,
		)
end


"""
`ir_mri_sensemap_sim_test1_show()`
"""
function ir_mri_sensemap_sim_test1_show(smap, x, y, zlist, title)
	clim = [-20,20]
	nz = length(zlist)
	pl = Array{Plot}(undef, nz)
	for iz = 1:nz
		pl[iz] = jim(x, y, smap[:,:,iz],
			title=title, clim=clim, xtick=[-2,2], ytick=[-2,2])
		blim = (zlist[iz] < 0.5) ? [7,12,19] : [1,3,5]
		contour!(pl[iz], x, y, abs.(smap[:,:,iz])', levels=blim,
			color=:blue, colorbar_entry=false)
		contour!(pl[iz], x, y, abs.(smap[:,:,iz])', levels=[0,0].+0.001,
			color=:green, colorbar_entry=false)
		title == "b" && plot!(xlabel="z = $(zlist[iz])")
	end
	return pl
end


"""
`ir_mri_sensemap_sim_test2(;chat)`
return plot with 2D example
"""
function ir_mri_sensemap_sim_test2( ; chat::Bool=true)

	@test_throws String ir_mri_sensemap_sim(scale=:bug)
	@test_throws String ir_mri_sensemap_sim(nring=2, orbit_start = [1,2,3])

	(smap,t) = ir_mri_sensemap_sim(dims=(32,32), scale=:ssos_center, chat=chat)

	if true # check rotational symmetry in 4-coil case
		for ic=2:4
			tmp = rotl90(smap[:,:,1], ic-1)
			@test isapprox(abs.(tmp), abs.(smap[:,:,ic]))
			p1 = angle.(tmp) .+ (ic-1) * pi/2 # add pi/2 to rotated
			p2 = angle.(smap[:,:,ic])
			tmp = ComplexF32.(cis.(p2 - p1))
			@test isapprox(tmp, ones(size(tmp))) # trick: equivs mod 2*pi
		end
	end

	ir_mri_sensemap_sim_show2(smap,
		t.x, t.y, t.dx, t.dy, t.nlist, t.plist, t.rlist)
end


"""
`ir_mri_sensemap_sim_test3(;chat)`
return plot that illustrates 3D sense maps
"""
function ir_mri_sensemap_sim_test3(;chat::Bool=false)
	nring = 3
	ncoil = 4 * nring
	ig = image_geom(nx=16, ny=14, nz=10, fov=200, dz=20, mask=:circ) # 20cm fov
#	ig = image_geom(nx=72, ny=48, nz=12, fov=22, zfov=10) % michelle

	(smap,t) = ir_mri_sensemap_sim(dims=(ig.nx, ig.ny, ig.nz),
			dx=ig.dx, dz=ig.dz,
			orbit_start = 1*[0,45,0],
			rcoil=70, nring=nring, ncoil=ncoil, chat=chat)

	if true
		tmp = smap .* repeat(ig.mask, 1,1,1,ncoil)
		jim(ncol=ig.nz, tmp, abswarn=false)
		chat && prompt()
		tmp = permutedims(tmp, [1,3,2,4]) # [nx nz ny ncoil] z cuts are smooth
		jim(ncol=ig.ny, tmp, abswarn=false)
		chat && prompt()
		jim(ncol=1, tmp[:,:,round(Int,end/2),:], abswarn=false)
		chat && prompt()
	end

	ir_mri_sensemap_sim_show3(smap, t.x, t.y, t.z, t.dx, t.dy, t.dz,
			t.nlist, t.plist, t.rlist, t.olist, t.ulist,
			t.nring, t.ncoilpr, t.rcoil)
end



"""
`ir_mri_sensemap_sim(:test)`
self test
"""
function ir_mri_sensemap_sim(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))

	@test typeof(ir_mri_sensemap_sim_test0()) <: Plots.Plot # ellipk
	@test typeof(ir_mri_sensemap_sim_test1()) <: Plots.Plot # basic test
	@test typeof(ir_mri_sensemap_sim_test2()) <: Plots.Plot # 2d test
	@test typeof(ir_mri_sensemap_sim_test3()) <: Plots.Plot # 3d test
	true
end

#= for plotting
ir_mri_sensemap_sim_test3(chat=true)
ir_mri_sensemap_sim_test3(chat=true)
=#
#@test ir_mri_sensemap_sim(:test)
