#=
disk-phantom.jl
generate random disk phantoms
2019-07-03, Jeff Fessler, University of Michigan
=#

export disk_phantom_params

# using MIRT: jim
using Random: seed!


"""
`params = disk_phantom_params( ; ...)`

Generate ellipse phantom parameters for a head-sized disk
plus many disks within it,
designed so that the disks have some minimum separation `minsep`
to avoid overlap and to simplify patch-based model fitting.
# Arguments
- `fov::Real = 240` image field of view in mm
- `rhead::Real = 100` background radius for "head" [mm]
- `muhead::Real = 1000` "mu" (intensity) value for background head disk
- `rmin::Real = 10` min radius for random disks
- `rmax::Real = 20` max radius for random disks
- `mumin::Real = 100` range of "mu" values for disks
- `mumax::Real = 300`
- `ndisk::Integer = 10` # of random disks
- `minsep::Real = 8` minimum disk separation in mm
- `maxtry::Integer = 500` give up on adding more disks if this is reached
- `warn::Bool = false` warn if maxtry reached?
- `seed::Integer = 0`if nonzero then use this seed
"""
function disk_phantom_params( ;
		fov::Real = 240,
		rhead::Real = 100, # background radius for "head" [mm]
		muhead::Real = 1000, # "mu" value for background head
		rmin::Real = 10, # min radius for disks
		rmax::Real = 20, # max radius for disks
		mumin::Real = 100, # range of "mu" values for disks
		mumax::Real = 300,
		ndisk::Integer = 10,
		minsep::Real = 8, # minimum separation in mm
		maxtry::Integer = 500, # give up on adding more disks if this is reached
		warn::Bool = false, # warn if maxtry reached?
		seed::Integer = 0,
	)

	params = zeros(Float32, ndisk+1, 6)
	params[end,:] = [0, 0, rhead, rhead, 0, muhead]

	randu = (a, b; f = x->x) -> a + f(rand()) * (b - a)
	seper = (a,b,r1,x,y,r2) -> sqrt.((x .- a).^2 + (y .- b).^2) .- r1 .- r2

	(seed != 0) && seed!(seed)

	for id=1:ndisk
		trial = 0
		for ii=1:maxtry
			rc = randu(0, rhead - rmin - minsep, f = x->x^0.5)
			phi = rand() * 2*pi
			rad = randu(rmin, rmax)
			val = randu(mumin, mumax)
			trial = [rc*cos(phi), rc*sin(phi), rad, rad, 0, val]

			id == 1 && break # first one is always fine

			# see if trial is too close to another one
			c = 1:(id-1) # check these
			sep = seper(params[c,1], params[c,2], params[c,3],
				trial[1], trial[2], trial[3])

			(minimum(sep) > minsep) && (rhead - (rc + rad) > minsep) && break
			ii == maxtry && warn && @warn("need more tries")
			ii == maxtry && return params[[1:id; end],:]
		end

#=
		if id > 1
			c = 1:(id-1) # check these
			sep = seper(params[c,1], params[c,2], params[c,3],
				trial[1], trial[2], trial[3])
			@show minimum(sep)
		end
		@show id, trial
=#

		params[id,:] = trial
	end

	return params
end


"""
`disk_phantom_params(:test)`
self test

`disk_phantom_params(:show)`
show example
"""
function disk_phantom_params(test::Symbol)
	ig = image_geom(nx = 128, fov = 240)
	params = disk_phantom_params( ; minsep = ig.dx*8)

	if test == :show
		tmp = ellipse_im(ig, params, oversample=3)
		return jim(ig.x, ig.y, tmp)
	end

	test != :test && throw(ArgumentError("test $test"))
	disk_phantom_params(:show)
	true
end

# disk_phantom_params(:test)
