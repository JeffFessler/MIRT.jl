# ellipse_im.jl

using MIRT: ellipse_im, ellipse_im_fast!, ellipse_im_params
using MIRT: image_geom, disk_phantom_params, shepp_logan_parameters
using MIRTjim: jim

using Plots: plot
#using Printf: @sprintf
using Test: @test, @test_throws, @inferred



"""
ellipse_im_show()
"""
function ellipse_im_show( ; over::Int = 2^2)
	ig = image_geom(nx=2^8, ny=2^8+2, fov=250)
#	ig.offset_y = 75.6 / ig.dy

	x0 = ellipse_im(ig, oversample=over)
	p0 = jim(ig.x, ig.y, x0, title="Shepp Logan", clim=(0.9,1.1))

	x1 = ellipse_im(ig, :shepplogan_emis, oversample=over)
	p1 = jim(ig.x, ig.y, x1, title="Shepp Logan Emission")

	x2 = ellipse_im(ig, :southpark)
	p2 = jim(ig.x, ig.y, x2, title="South Park")

	x3 = ellipse_im(ig, :shepplogan_brainweb)
	p3 = jim(ig.x, ig.y, x3, title="Shepp Logan Brainweb")

	plot(p0,p1,p2,p3)
end


#=

# compare to aspire
function ellipse_im_aspire()
	nx = 2^6
	ig = image_geom(nx=nx, ny=nx+2, fov=2^7)
	ell = [10, 20, 30, 40, 50, 1]

	file = tempname() * "t.fld"
	file = "tt.fld"
	over = 2^2
#	com = @sprintf "echo y | op -chat 99 ellipse %s %d %d  %g %g %g %g %g %g %d" file ig.nx ig.ny ell ./ [ig.dx, ig.dx, ig.dx, ig.dx, 1, 1] log2(over)+1
	tmp = ""
	for i=1:4
		tmp = tmp * " $(ell[i] / ig.dx)"
	end
	for i=5:6
		tmp = tmp * " $(ell[i])"
	end
	op = "/Users/fessler/bin/mi64-gcc/op"
	com = "$op -chat 99 ellipse $file $(ig.nx) $(ig.ny)$tmp $(log2(over)+1)"
	com = "/bin/ls /Users/fessler/bin/mi64-gcc" # todo: why fails?
	@show com
#	cm = @cmd "$com"
	cm = `"$com"`
	@show cm

	run(cm)
return true
	asp = fld_read(file)
	p1 = jim(asp, title="aspire")

	area_asp = sum(asp) * abs(ig.dx * ig.dy)

	jul = ellipse_im(ig, ell, oversample=over) # 'type', types{ii})
	area_jul = sum(jul) * abs(ig.dx * ig.dy)

	p2 = jim(jul, title="jul")
	p3 = jim((jul-asp)*over^2, "difference: (jul-asp)*over^2")

	@show maximum(abs.(jul - mat)) / ell[6] * over^2
	!isapprox(jul, mat) && throw("approximation error")
#	max_percent_diff(mat, asp)

	area_real = pi * ell[3] * ell[4] * ell[6]
	@show area_real, area_asp, area_jul

	true
end
=#


#=
% ellipse_im_profile()
function ellipse_im_profile
ig = image_geom('nx', 2^9, 'ny', 2^9+2', 'fov', 250)
profile on
x0 = ellipse_im(ig, [], 'oversample', 3, 'type', 'fast')
profile off
profile report
=#


#@test ellipse_im_aspire() # todo

fov = 100
@inferred shepp_logan_parameters(fov, fov)
shepp_logan_parameters(fov, fov, case=:kak)
shepp_logan_parameters(fov, fov, case=:emis)
shepp_logan_parameters(fov, fov, case=:brainweb)

# test various ways of calling
ellipse_im(20)
#@inferred ellipse_im(20,22) # todo-i: why fails?
ellipse_im(20,22)
ellipse_im(30, :shepplogan_emis, oversample=2)

ig = image_geom(nx=80, dx=1)
ellipse_im_params(ig, :kak)

params = shepp_logan_parameters(ig.fovs..., case=:brainweb)
phantom = zeros(Float32, ig.dims)
@inferred ellipse_im_fast!(phantom, ig.nx, ig.ny, params, ig.dx, ig.dy,
		ig.offset_x, ig.offset_y, 0, 2, false, 1.)
ellipse_im(ig, :disks)
ellipse_im(ig, params, oversample=2)
ellipse_im(ig, params, how=:fast, replace=true, rot=30)

@test_throws String ellipse_im(ig, :bad)
@test_throws String ellipse_im(ig, params, how=:bad)
@test_throws String ellipse_im(ig, params')
@test_throws String shepp_logan_parameters(ig.fovs..., case=:bad)

ellipse_im(ig, params)

#=
x1 = ellipse_im(100, :shepplogan_emis)
x2 = ellipse_im(100, :shepplogan_emis, oversample=2)
plot(jim(x1), jim(x2), jim(x2-x1))
=#

ellipse_im_show()
