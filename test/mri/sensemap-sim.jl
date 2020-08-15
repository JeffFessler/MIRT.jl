# sensemap-sim.jl

using MIRT: ir_mri_sensemap_sim
using MIRT: jim, prompt, image_geom, ndgrid
import MIRT: ir_mri_smap_r
import MIRT: ir_mri_smap1

using SpecialFunctions: ellipk, ellipe
using Test: @test, @testset, @test_throws, @inferred
using Plots: Plot, plot, plot!, contour!, scatter!, quiver, quiver!, Arrow


# many plots
function sensemap_test()
	nx,ny = 32,32
	smap, info = ir_mri_sensemap_sim(:all ;
		dims=(nx,ny), scale=:ssos_center,
	)
	smap = reshape(smap, nx, ny, 1, info.ncoilpr, info.nring)

	x = info.x
	y = info.y
	data = info.data
	down2 = a -> vec(a[1:3:end,1:3:end])
	down3 = a -> vec(a[1:3:end,1:3:end,1:3:end])

	for ir = 1:info.nring
		for ic = 1:info.ncoilpr
			pc = data[ir,ic]
			xr = pc.xr
			zr = pc.zr

			p = Array{Any}(undef, 7)

			# coordinates
			p[1] = jim(x, y, xr, xlabel="x", ylabel="y", "xr")
			p[2] = jim(x, y, zr, "zr")

			sx = pc.sx
			sy = pc.sy
			sz = pc.sz
			tmp = sqrt.(sx.^2 + sz.^2)
			(xx,yy) = ndgrid(x,y)

			p[3] = jim(x, y, sx, "sx")
			p[4] = jim(x, y, sy, "sy")
			p[5] = jim(x, y, sz, "sz")
			p[6] = quiver(down2(xx), down2(yy), title="(sx,sy)",
					quiver=(down2(sx./tmp), down2(sz./tmp)))

			# see final field components vs phase (for nz=0)
			(xx,yy) = ndgrid(x,y)
			bx = pc.bx
			by = pc.by
			bb = sqrt.(bx.^2 + by.^2)
			p[7] = jim(x, y, angle.(smap[:,:,1,ic,1]), "phase", color=:hsv)
			quiver!(down2(xx), down2(yy), title="(bx,by)",
					quiver=(down3((bx./bb)), down3((by./bb))))
			plot(p...)
			prompt()
		end
	end

end


# test code to explore when r is near 0
function ir_mri_smap_r_test()
	r0 = LinRange(0,5e-7,101)
	z0 = 0.4
	t0 = ir_mri_smap_r.(r0, z0)
	slope = 3*pi * z0 / ((1+z0^2)^2.5)
	plot(r0, t0, label="fun", xlabel="r")
	plot!(r0, slope * r0, label="line")
end


"""
ellipke_plot()
"""
function ellipke_plot()
	m = LinRange(0,1,101)
	(k,e) = (ellipk.(m), ellipe.(m))
	plot(xaxis=[0,1], yaxis=[0,3π/2])
	plot!(m, k, label="k", xlabel="m")
	plot!(m, e, label="e")
	plot!(ytick=((0:3)*pi/2, ["0", "π/2", "π", "3π/2"]))
end


"""
    ir_mri_sensemap_sim_test1()
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
    ir_mri_sensemap_sim_test1_show()
"""
function ir_mri_sensemap_sim_test1_show(smap, x, y, zlist, title)
	clim = (-20,20)
	nz = length(zlist)
	pl = Array{Any}(undef, nz)
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
    ir_mri_sensemap_sim_show2()
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

	pl = Array{Any}(undef, ncoil, 3)
	clim = (0, maximum(abs.(smap)))
	xmax = maximum(abs.([vec(x); vec(y); vec(plist[:,:,[1,2]])]))
	ymax = xmax
	for ic=1:ncoil
		tmp = smap[:,:,ic]
		p = jim(x, y, abs.(tmp), clim=clim, "Magn. $ic")
		plot!(p, xlim=[-1,1]*1.1*xmax, xtick=(-1:1) * nx/2 * dx)
		plot!(p, ylim=[-1,1]*1.1*xmax, ytick=(-1:1) * nx/2 * dy)

		scatter!(p, [0], [0], marker=:o, label="", color=:green) # center
		scatter!(p, vec(plist[:,:,1]), vec(plist[:,:,2]),
			marker=:o, label="", color=:blue) # coil location
		xdir = nlist[ic,1,2]
		ydir = nlist[ic,1,1]
		r = rlist[ic,1]
		plot!(p, plist[ic,1,1].+r*xdir*[-1,1], plist[ic,1,2].+r*ydir*[1,-1],
			label="", color=:blue, linewidth=3) # coil
		pl[ic,1] = p

		ph = angle.(tmp) # show raw phase (understandable with hsv colormap)
		p = jim(x, y, ph, "Phase", clim=(-pi,pi), color=:hsv)
		plot!(p, xlim=[-1,1]*1.1*xmax, xtick=(-1:1) * nx/2 * dx)
		plot!(p, ylim=[-1,1]*1.1*xmax, ytick=(-1:1) * nx/2 * dy)
		pl[ic,2] = p
	end

	ssos = sqrt.(sum(abs2.(smap), dims=ndims(smap)))
	ssos = ssos / ssos[Int(end/2),Int(end/2)]

	p = jim(x, y, ssos, "SSoS (norm.)",
		#	xlim=[-1,1]*1.1*xmax,
			xtick=(-1:1) * nx/2 * dx,
			ytick=(-1:1) * ny/2 * dy)
	pl[1,3] = p

	if true # quiver plot for 1st coil
		bx = real(smap[:,:,1])
		by = imag(smap[:,:,1])
		(xx,yy) = ndgrid(x,y)
		down = a -> vec(a[1:3:end,1:3:end])
		p = quiver(down(xx), down(yy), quiver=(down(bx), down(by)),
				aspect_ratio = 1,
				arrow = Arrow(:simple, :head, 0.01, 0.01), # no effect!?
				title = "Field in x-y plane")
		plot!(p, xlim=[-1,1]*1.1*xmax, xtick=(-1:1) * nx/2 * dx)
		plot!(p, ylim=[-1,1]*1.1*xmax, ytick=(-1:1) * nx/2 * dy)
		pl[2,3] = p
	end

	for ic=3:ncoil
		pl[ic,3] = plot(xaxis=:off,yaxis=:off,grid=:off) # kludge
	end

	return plot(pl...)
end


"""
    ir_mri_sensemap_sim_test2( ; chat)
return plot with 2D example
"""
function ir_mri_sensemap_sim_test2()

	@test_throws String ir_mri_sensemap_sim(scale=:bug)
	@test_throws String ir_mri_sensemap_sim(nring=2, orbit_start = [1,2,3])

	(smap,t) = ir_mri_sensemap_sim(
		Vector{Tuple{Int,Int}}(undef, 0) ;
		dims=(32,32), scale=:ssos_center,
	)

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

	return ir_mri_sensemap_sim_show2(smap,
		t.x, t.y, t.dx, t.dy, t.nlist, t.plist, t.rlist)
end


"""
    ir_plot3_cube(x,y,z)
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
    ir_mri_sensemap_sim_show3()
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
    ir_mri_sensemap_sim_test3( ; chat)
return plot that illustrates 3D sense maps
"""
function ir_mri_sensemap_sim_test3( ; chat::Bool=false)
	nring = 3
	ncoil = 4 * nring
	ig = image_geom(nx=16, ny=14, nz=10, fov=200, dz=20, mask=:circ) # 20cm fov
#	ig = image_geom(nx=72, ny=48, nz=12, fov=22, zfov=10) % michelle

	(smap,t) = #@inferred # todo-i fails
		ir_mri_sensemap_sim(
			Vector{Tuple{Int,Int}}(undef, 0) ;
			dims=(ig.nx, ig.ny, ig.nz),
			dx=ig.dx, dz=ig.dz,
			orbit_start = 1*[0,45,0],
			rcoil=70, nring=nring, ncoil=ncoil,
		)

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


@testset "1D" begin
	@inferred ir_mri_smap_r(5e-7, 0.4)
	plot(
		ir_mri_smap_r_test(),
		ellipke_plot(),
	)
end
prompt()

@testset "smap1" begin
	@inferred ir_mri_smap1(1, 2, 3, 4)
	ir_mri_sensemap_sim_test1() # basic test
end
prompt()

#@inferred ir_mri_sensemap_sim() # fails: Union of 3D, 4D

ir_mri_sensemap_sim_test2() # 2d test
prompt()

ir_mri_sensemap_sim_test3() # 3d test
prompt()

sensemap_test()

#=
@test typeof(ir_mri_sensemap_sim_test0()) <: Plots.Plot # ellipk
@test typeof(ir_mri_sensemap_sim_test1()) <: Plots.Plot # basic test
@test typeof(ir_mri_sensemap_sim_test2()) <: Plots.Plot # 2d test
@test typeof(ir_mri_sensemap_sim_test3()) <: Plots.Plot # 3d test
=#

#= for plotting
ir_mri_sensemap_sim_test3(chat=true)
ir_mri_sensemap_sim_test3(chat=true)
=#
