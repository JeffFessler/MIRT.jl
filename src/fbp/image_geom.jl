#=
image_geom.jl
Methods related to an image geometry for image reconstruction
2017-10-02 Samuel Rohrer, University of Michigan
2019-03-05 Jeff Fessler, Julia 1.1 + tests
2019-06-23 Jeff Fessler, overhaul
=#

export MIRT_image_geom, image_geom, cbct

#using MIRT: jim, downsample2, downsample3
using Test: @test
using ImageTransformations: imresize


"""
`MIRT_image_geom`

Image geometry "struct" with essential parameters
"""
struct MIRT_image_geom
	# options for 2D image geometry
	nx::Int				# image dimension 1
	ny::Int				# image dimension 2
	dx::Float32			# pixel size
	dy::Float32			# pixel size
	offset_x::Float32	# unitless
	offset_y::Float32	# unitless

	# options for 3D image geometry
	nz::Int				# image dimension 3
	dz::Float32			# voxel size
	offset_z::Float32	# unitless

	mask::Array{Bool}	# logical mask
end


"""
`image_geom_help( ; io)`
"""
function image_geom_help( ; io::IO = isinteractive() ? stdout : IOBuffer())
	print(io, "propertynames:\n")
	print(io, propertynames(image_geom(nx=1, dx=1)))

	print(io,
	"
	Derived values for 2D (and 3D) images:

	is3	is it 3d (nz > 0?)
	fovs [|dx|*nx |dy|*ny ...]
	np	sum(mask) = # pixels to be estimated

	dim	(nx, ny, [nz])
	x	1D x coordinates of each pixel
	y	1D y coordiantes of each pixel
	wx	(nx - 1/2) * dx + offset_x
	wy	(ny - 1/2) * dy + offset_y
	wz	(nz - 1/2) * dz + offset_z

	xg	x coordinates of each pixel as a grid (2D or 3D)
	yg	y coordinates of each pixel as a grid (2D or 3D)

	u	1D frequency domain coordinates of each pixel (kx)
	v	1D frequency domain coordinates of each pixel (ky)
	w	1D frequency domain coordinates of each pixel (kz)

	ug	2D or 3D grid of frequency domain coordinates
	vg	2D or 3D grid of frequency domain coordinates
	fg	NamedTuple of 2D or 3D frequency coordinates
	ones	ones(dim)
	zeros	zeros(dim)

	mask_outline	binary image showing 2D outline of mask

	Derived values for 3D images:

	z	z coordinates of each pixel
	zg	z coordinates of each pixel as a grid (3D only)
	wg	3D grid of frequency domain coordinates

	mask_or	[nx ny] logical 'or' of mask

	Methods:

	embed(x)	turn short column(s) into array(s)
	maskit(x)	opposite of embed
	shape(x)	reshape long column x to [nx ny [nz]] array
	unitv(...)	j=j | i=(ix,iy) | c=[cx cy]
				j: single index from 1 to length(z)
				i: [ix,iy[,iz]] index from 1 to nx,ny
				c: [cx,cy[,cz]] index from +/- n/2 center at floor(n/2)+1
	circ(rx=,ry=,cx=,cy=)	circle of given radius and center (cylinder in 3D)
	plot()		plot the image geometry

	Methods that return a new MIRT_image_geom:

	down(down::Int)		down-sample geometry by given factor
	over(over::Int)		over-sample geometry by given factor
	expand_nz(nz_pad)	expand image geometry in z by nz_pad on both ends
	\n")
end


"""
`MIRT_cbct_ig`

Structure suitable for passing to C routines `cbct_*`
based on the struct `cbct_ig` found in `cbct,def.h`
"""
struct MIRT_cbct_ig
	nx::Cint
	ny::Cint
	nz::Cint
	dx::Cfloat
	dy::Cfloat
	dz::Cfloat
	offset_x::Cfloat
	offset_y::Cfloat
	offset_z::Cfloat
	mask2::Ptr{Cuchar}	# [nx,ny] 2D support mask: 0 or 1 .. nthread
	iy_start::Ptr{Cint} # [nthread] for mask2
	iy_end::Ptr{Cint}	# [nthread] for mask2
end

"""
`cbct(ig::MIRT_image_geom; nthread=1)`
constructor for `MIRT_cbct_ig`
"""
function cbct(ig::MIRT_image_geom; nthread::Int=1)
	iy_start = [0]
	iy_end = [ig.ny]
	nthread != 1 && throw("only nthread=1 for now due to iy_start")
	return MIRT_cbct_ig(
		Cint(ig.nx), Cint(ig.ny), Cint(ig.nz), Cfloat(ig.dx),
		Cfloat(ig.dy), Cfloat(ig.dz),
		Cfloat(ig.offset_x), Cfloat(ig.offset_y), Cfloat(ig.offset_z),
			pointer(UInt8.(ig.mask_or)),
			pointer(Cint.(iy_start)), pointer(Cint.(iy_end)),
		)
end


# see https://docs.julialang.org/en/stable/manual/constructors/
"""
`ig = image_geom(...)`

Constructor for `MIRT_image_geom`

# Arguments
- `nx::Int = 128`
- `ny::Int = nx`
- `dx::Real = ?` (must specify one of `dx` or `fov`)
- `dy::Real = -dx`
- `offset_x::Real = 0` (unitless)
- `offset_y::Real = 0` (unitless)
- `fov::Real = ?` (if specified, then `nx*dx=ny*dy`)
- `nz::Int = 0`
- `dz::Real = ?` (need one of `dz` or `zfov` if `nz > 0`)
- `zfov::Real = ?` (if specified, then `nz*dz`)
- `offset_z::Real = 0` (unitless)
- `offsets::Symbol = :none` or :dsp
- `mask::Union{Symbol,AbstractArray{Bool}} = :all` | `:circ` | `:all_but_edge`
"""
function image_geom( ;
		nx::Integer = 128,
		ny::Integer = nx,
		dx::Real = NaN,
		dy::Real = NaN,
		offset_x::Real = 0,
		offset_y::Real = 0,
		fov::Real = NaN,
		nz::Integer = 0,
		dz::Real = NaN,
		zfov::Real = NaN,
		offset_z::Real = 0,
		offsets::Symbol = :none,
		mask::Union{Symbol,AbstractArray{Bool}}	= :all,
	)

	# handle optional arguments

	is3 = !isnan(dz) || !isnan(zfov) || (nz > 0)

	# offsets
	if offsets == :dsp
		# check that no nonzero offset was specified
		(offset_x != 0 || offset_y != 0) && throw("offsets usage incorrect")
		offset_x = 0.5
		offset_y = 0.5
		is3 && offset_z != 0 && throw("offsets usage incorrect")
	elseif offsets != :none
		throw("offsets $offsets")
	end

	# transverse plane (2d) distances
	if true
		if isnan(fov)
			isnan(dx) && throw("dx or fov required")
		#	fov = nx * dx
		else
			!isnan(dx) && throw("dx and fov?")
			dx = fov / nx
			dy = fov / ny
		end

		dy = isnan(dy) ? -dx : dy # default dy
	end

	# 3D geometry
	if is3
		if isnan(zfov)
			isnan(dz) && throw("dz or zfov required")
		#	zfov = nz * dz
		else
			!isnan(dz) && throw("dz and zfov?")
			dz = zfov / nz
		end
	end

	# mask
	if mask == :all
		mask = is3 ? trues(nx, ny, nz) : trues(nx, ny)
	elseif mask == :circ
		mask = image_geom_circle(nx, ny, dx, dy)
	elseif mask == :all_but_edge
		mask = trues(nx,ny,nz)
		mask[  1,   :, :] .= false
		mask[end,   :, :] .= false
		mask[  :,   1, :] .= false
		mask[  :, end, :] .= false
	elseif isa(mask, Symbol)
		throw("mask symbol $mask")
	elseif size(mask,1) != nx || size(mask,2) != ny || (is3 && size(mask,3) != nz)
		throw("mask size $(size(mask)), nx=$nx ny=$ny nz=$nz")
	end

	# return the object
	return MIRT_image_geom(nx, ny, dx, dy, offset_x, offset_y,
		nz, dz, offset_z, mask)
end


"""
`image_geom_down_round()`
helper function needed to downsample `image_geom`
"""
function _down_round(val::Real=1, dd::Real=1, down::Real=1)
	out = val / down
	if out != round(out)
		out = 2 * round(out / 2) # keep it even if not an integer
	end
	dd = dd * down
	return Int(out), dd
end


"""
`ig_down = downsample(ig, down::Union{Int,Vector{Int}})`
down sample an image geometry by the factor `down`
cf `image_geom_downsample`
"""
function downsample(ig::MIRT_image_geom, down::Union{Int,Vector{Int}})

	# check if it is a scalar by seeing if the collection of size is empty
	if isempty(size(down))
		downv = ig.is3 ? [down,down,down] : downv = [down,down]
	else
		downv = down
	end

	# call the down round function
	down_nx, down_dx = _down_round(ig.nx, ig.dx, downv[1])
	down_ny, down_dy = _down_round(ig.ny, ig.dy, downv[2])
	# adjust to "pixel" units
	down_offset_x = ig.offset_x / downv[1]
	down_offset_y = ig.offset_y / downv[2]

	if ig.is3 # 3d case
		down_nz, down_dz = _down_round(ig.nz, ig.dz, downv[3])
		down_offset_z = ig.offset_z / downv[3]
	else
		down_nz = 0; down_dz = 0; down_offset_z = 0
	end

	# carefully down-sample the mask
	mdim_vec = [size(ig.mask)...] # tuple to vector
	if ig.is3
		if [down_nx, down_ny, down_nz] .* downv == mdim_vec
			down_mask = downsample3(ig.mask, downv) .> 0
		else
			throw("bug: bad mask size. need to address mask downsampling")
		end
	else
		if [down_nx,down_ny] .* downv == mdim_vec
			down_mask = downsample2(ig.mask, downv) .> 0
		else
			down_mask = imresize(ig.mask, (down_nx,down_ny)...) .> 0
		end
	end

	return MIRT_image_geom(
		down_nx, down_ny, down_dx, down_dy, down_offset_x, down_offset_y,
		down_nz, down_dz, down_offset_z, down_mask)
end


"""
`ig_new = image_geom_expand_nz(ig::MIRT_image_geom, nz_pad::Integer)`
pad both ends
"""
function image_geom_expand_nz(ig::MIRT_image_geom, nz_pad::Integer)
	!ig.is3 && throw("expand_nz only valid for 3D")
	out_nz = ig.nz + 2*nz_pad
	out_mask = cat(dims=3, repeat(ig.mask[:,:,1], 1, 1, nz_pad), ig.mask,
			repeat(ig.mask[:,:,end], 1, 1, nz_pad))
	return MIRT_image_geom(
		ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y,
		out_nz, ig.dz, ig.offset_z, out_mask)
end


"""
`ig_over = image_geom_over(ig::MIRT_image_geom, over::Integer)`
over-sample an image geometry by the factor `over`
"""
function image_geom_over(ig::MIRT_image_geom, over::Integer)
	if all(ig.mask .== true)
		mask_over = trues(ig.dim)
	else
		mask_over = imresize(ig.mask, ig.dim .* over) .> 0
	end
	return MIRT_image_geom(
		ig.nx*over, ig.ny*over, ig.dx/over, ig.dy/over,
		ig.offset_x*over, ig.offset_y*over,
		ig.nz*over, ig.dz/over, ig.offset_z*over,
		mask_over)
end


# default is a circle that just inscribes the square
# but keeping a 1 pixel border due to ASPIRE regularization restriction
function image_geom_circle(nx::Int, ny::Int, dx::Real, dy::Real;
	rx::Real = min(abs((nx/2-1)*dx), abs((ny/2-1)*dy)),
	ry::Real = rx, cx::Real = 0, cy::Real = 0, nz::Int=0, over::Int=2)
	ig = image_geom(nx=nx, ny=ny, dx=dx, dy=dy)
	circ = ellipse_im(ig, [cx cy rx ry 0 1], oversample=over) .> 0
	if nz > 0
		circ = repeat(circ, 1, 1, nz)
	end
	return circ
end


"""
`out = image_geom_add_unitv(z; j=?, i=?, c=?)`

add a unit vector to an initial array `z` (typically of zeros)

# Arguments
- `j` single index from 1 to length(z)
- `i` [ix,iy[,iz]] index from 1 to nx,ny
- `c` [cx,cy[,cz]] index from +/- n/2 center at floor(n/2)+1

default with no arguments gives unit vector at center `c=[0,0]`
"""
function image_geom_add_unitv(z; # starts with zeros()
		j::Integer=0,
		i::AbstractVector{<:Int} = zeros(Int, ndims(z)),
		c::AbstractVector{<:Int} = zeros(Int, ndims(z))
		)
	out = copy(z)

	if j > 0 && all(i .== 0)
		out[j] += 1
	elseif j == 0 && any(i .> 0)
		out[i...] += 1
	else
		tmp = c .+ Int.(floor.(size(out) ./ 2)) .+ 1
		out[tmp...] += 1
	end

	return out
end


"""
`image_geom_add_unitv(:test)`
"""
function image_geom_add_unitv(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	image_geom_add_unitv(zeros(3,4), j=2)
	image_geom_add_unitv(zeros(3,4), i=[2,3])
	image_geom_add_unitv(zeros(3,4))
	image_geom_add_unitv(zeros(3,4), c=[-1,-2])
	true
end


"""
`image_geom_plot(ig)`
"""
function image_geom_plot(ig::MIRT_image_geom; kwargs...)
	return ig.is3 ?
		jim(ig.x, ig.y, ig.mask_or,
			"(dx,dy,dz)=$(ig.dx),$(ig.dy),$(ig.dz)"; kwargs...) :
		jim(ig.x, ig.y, ig.mask, "(nx,ny)=$(ig.nx),$(ig.ny)"; kwargs...)
end


function Base.display(ig::MIRT_image_geom)
	ir_dump(ig)
end


# Extended properties

image_geom_fun0 = Dict([
	(:help, ig -> image_geom_help()),

	(:is3, ig -> ig.nz > 0),
	(:dim, ig -> ig.is3 ? (ig.nx, ig.ny, ig.nz) : (ig.nx, ig.ny)),
	(:fovs, ig -> ig.is3 ?
		[abs(ig.dx)*ig.nx, abs(ig.dy)*ig.ny, abs(ig.dz)*ig.nz] :
			[abs(ig.dx)*ig.nx, abs(ig.dy)*ig.ny]),

	(:zeros, ig -> zeros(Float32, ig.dim)), # (nx, ny, nz) : zeros(nx, ny)
	(:ones, ig -> ones(Float32, ig.dim)), # (nx, ny, nz) : zeros(nx, ny)

	(:x, ig -> ((0:ig.nx-1) .- ig.wx) * ig.dx),
	(:y, ig -> ((0:ig.ny-1) .- ig.wy) * ig.dy),
	(:z, ig -> ((0:ig.nz-1) .- ig.wz) * ig.dz),

	(:wx, ig -> (ig.nx - 1)/2 + ig.offset_x),
	(:wy, ig -> (ig.ny - 1)/2 + ig.offset_y),
	(:wz, ig -> (ig.nz - 1)/2 + ig.offset_z),

	(:xg, ig -> ig.is3 ? repeat(ig.x, 1, ig.ny, ig.nz) :
						repeat(ig.x, 1, ig.ny)),
	(:yg, ig -> ig.is3 ? repeat(ig.y', ig.nx, 1, ig.nz) :
						repeat(ig.y', ig.nx, 1)),
	(:zg, ig -> ig.is3 ? repeat(reshape(ig.z, 1, 1, ig.nz), ig.nx, ig.ny, 1) :
						zeros(ig.nx, ig.ny)),

	# DFT frequency sample grid
	(:u, ig -> (-ig.nx/2:(ig.nx/2-1)) / (ig.nx*ig.dx)),
	(:v, ig -> (-ig.ny/2:(ig.ny/2-1)) / (ig.ny*ig.dy)),
	(:w, ig -> (-ig.nz/2:(ig.nz/2-1)) / (ig.nz*ig.dz)),

	(:ug, ig -> ig.is3 ? repeat(ig.u, 1, ig.ny, ig.nz) :
						repeat(ig.u, 1, ig.ny)),
	(:vg, ig -> ig.is3 ? repeat(ig.v', ig.nx, 1, ig.nz) :
						repeat(ig.v', ig.nx, 1)),
	(:wg, ig -> ig.is3 ? repeat(reshape(ig.w, 1, 1, ig.nz), ig.nx, ig.ny, 1) :
						zeros(ig.nx, ig.ny)),
	(:fg, ig -> ig.is3 ? (ig.ug, ig.vg, ig.wg) : (ig.ug, ig.vg)),

	(:np, ig -> sum(ig.mask)),
	(:mask_or, ig -> mask_or(ig.mask)),
	(:mask_outline, ig -> mask_outline(ig.mask)),

	# simple functions

	(:plot, ig -> ((;kwargs...) -> image_geom_plot(ig; kwargs...))),
	(:embed, ig -> (x::AbstractArray{<:Number} -> embed(x, ig.mask))),
	(:maskit, ig -> (x::AbstractArray{<:Number} -> maskit(x, ig.mask))),
	(:shape, ig -> (x::AbstractArray{<:Number} -> reshape(x, ig.dim))),
	(:unitv, ig -> ((;kwargs...) -> image_geom_add_unitv(ig.zeros; kwargs...))),
	(:circ, ig -> ((;kwargs...) ->
		image_geom_circle(ig.nx,ig.ny,ig.dx,ig.dy,nz=ig.nz; kwargs...))),

	# functions that return new geometry:

	(:down, ig -> (down::Int -> downsample(ig, down))),
	(:over, ig -> (over::Int -> image_geom_over(ig, over))),
	(:expand_nz, ig -> (nz_pad::Int -> image_geom_expand_nz(ig, nz_pad))),

	])


# Tricky overloading here!

Base.getproperty(ig::MIRT_image_geom, s::Symbol) =
		haskey(image_geom_fun0, s) ? image_geom_fun0[s](ig) :
		getfield(ig, s)

Base.propertynames(ig::MIRT_image_geom) =
	(fieldnames(typeof(ig))..., keys(image_geom_fun0)...)


function image_geom_test2(ig::MIRT_image_geom)
	# test 2D functions provided by the constructor
	ig.dim
	ig.x
	ig.y
	ig.wx
	ig.wy
	ig.xg
	ig.yg
	ig.fovs
	ig.np
	ig.mask_outline
	ig.ones
	ig.zeros
	ig.u
	ig.v
	ig.ug
	ig.vg
	ig.fg
	@test ig.shape(ig.ones[:]) == ig.ones
	@test ig.embed(ig.ones[ig.mask]) == Float32.(ig.mask)
	@test ig.maskit(ig.ones) == Float32.(ones(ig.np))
	ig.unitv(j=4)
	ig.unitv(i=ones(Int, length(ig.dim)))
	ig.unitv(c=zeros(Int, length(ig.dim)))
	ig.unitv()
	ig.circ()
	ig.plot()
	ig.down(2)
	ig.over(2)
	true
end

function image_geom_test2()
	image_geom(nx=16, dx=2, offsets=:dsp, mask=:all_but_edge)
	@test_throws String image_geom(nx=16, dx=1, offsets=:bad)
	@test_throws String image_geom(nx=16, dx=1, mask=:bad)
	ig = image_geom(nx=16, dx=2)
	display(ig)
	image_geom_test2(ig)
	ig = image_geom(nx=16, dx=2, mask=:circ)
	ig.over(2)
	ig.down(3) # test both even and non-even factors
	ig.help
	true
end


function image_geom_test3(ig::MIRT_image_geom)
	@test image_geom_test2(ig)
	ig.wz
	ig.zg
	ig.mask_or
	ig.expand_nz(2)
	cbct(ig)
	true
end

function image_geom_test3()
	ig = image_geom(nx=16, nz=4, dx=2, zfov=1)
	ig = image_geom(nx=16, nz=4, dx=2, dz=3)
	image_geom_test3(ig)
	true
end


"""
`image_geom(:test)`
self test
"""
function image_geom(test::Symbol)
	if test == :help
		image_geom_help()
		return true
	end
	test != :test && throw(ArgumentError("test $test"))
	@test _down_round(4, 3, 2)[1] == 2
	@test image_geom_test2()
	@test image_geom_test3()
	@test image_geom(:help)
	true
end
