# image_geom.jl
# Methods related to an image geometry for image reconstruction
# 2017-10-02 Samuel Rohrer, University of Michigan
# 2019-03-05 Jeff Fessler, Julia 1.1 + tests

#using Printf
using Test
using DSP: conv2
#using SparseArrays: sparse, findnz, AbstractSparseVector


"""
`array = embed(v, mask)`

embed vector `v` of length `sum(mask)`
into elements of an array where `mask` is `true`
"""
function embed(
		v::AbstractVector{<:Number},
		mask::AbstractArray{Bool,N} where N)
	array = zeros(eltype(v), size(mask))
	array[mask] .= v
	return array
end


#=
this is needed only for "unpacking" sparse system matrices so ignore for now
# image_geom_embed_sparse
# called by image_geom_embed
# function _embed_sparse(x::Array{T}) where {T <: Number}
function embed(x::AbstractSparseVector{<:Number},
				mask::AbstractArray{Bool,N} where N)
	i, v = findnz(x)
	ind = findall(mask)
	j = ind(j)
	return sparsevec(i, j, a, size(x,1), length(mask))
end
=#


"""
`embed(:test)`
"""
function embed(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask = [false true true; true false false]
	@test embed(1:3,mask) == [0 2 3; 1 0 0]
#	@test embed(sparse(1:3),mask) == sparse([0 2 3; 1 0 0]) # later
	true
end


"""
`maskit(x::AbstractArray{<:Number})`
opposite of embed
"""
function maskit(x::AbstractArray{<:Number}, mask::Array{Bool})
	dim = size(x)
	x = reshape(x, length(mask), :)
	x = x[mask[:],:] # reshape(mask, prod(_dim()))]
	if length(dim) == ndims(mask)
		x = dropdims(x, dims=2) # squeeze
	elseif length(dim) > ndims(mask)
		x = reshape(x, :, dim[(1+ndims(mask)):end])
	else
		throw(DimensionMismatch("size(x) = $(size(x))"))
	end
	return x
end


"""
`out = image_geom_add_unitv(z; j=?, i=?, c=?)`

add a unit vector to an initial array `z` (typically of zeros)

option;
`j` single index from 1 to length(z)
`i` [ix,iy[,iz]] index from 1 to nx,ny
`c` [cx,cy[,cz]] index from +/- n/2 center at floor(n/2)+1

default with no arguments gives unit vector at center c=[0,0]
"""
function image_geom_add_unitv(z; # starts with zeros()
		j::Integer=0,
		i::AbstractVector{<:Int} = zeros(Int, ndims(z)),
		c::AbstractVector{<:Int} = zeros(Int, ndims(z))
		) # todo where {T1 <: Integer, T2 <: Integer}
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


# Image geometry "struct" with some parameters and methods included
struct MIRT_image_geom
  # options for 2D image geometry
  nx::Integer		# image dimension one
  ny::Integer		# image dimension (default: nx) two
  dx::Real			# pixel size ('dx' or 'fov' required)
  dy::Real			# pixel size (default: -dx). (value | 'dx' | '-dx' | '2*dx')
  offset_x::Real	# unitless (default: 0)
  offset_y::Real	# unitless (default: 0)
  fov::Real			# nx*dx = ny*dy (if specified)

  # options for 3D image geometry
  nz::Integer		# image dimension three
  dz::Real			# pixel size (default: dx)
  zfov::Real		# dz*nz
  offset_z::Real	# unitless (default: 0)

  # general image orientation
  offsets::String	# this is a string, only 'dsp' has meaning now
  mask_type::String	# a string stored to make downsample easier
  mask::Array{Bool}	# logical mask
  iy_start::Array{Int}	# size [nthread]
  iy_end::Array{Int}	# size [nthread]
  is3::Bool			# storing this to make life easier

  # functions provided by the constructor
  # methods for 2D images
#  unitv::Function		# (ix,iy) | (jj) | ('c', [cx cy])
  circ::Function		# (rx,ry,cx,cy) circle of radius rx, ry (cylinder in 3D)
  expand_nz::Function	# '(nx_pad)'
  over::Function		#
  plot::Function		# plot the image geometry

  # methods for 3D images
end


"""
`image_geom_help()`
"""
function image_geom_help()
	print("propertynames:\n")
	print(propertynames(image_geom(nx=1,dx=1)))

	"
	Derived values for 2D (and 3D) images:

	fovs [|dx|*nx |dy|*ny ...]
	np	sum(mask) = # pixels to be estimated

	dim	[nx ny [nz]]
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
	"
end


# structure suitable for passing to C routines cbct_*
# based on the struct cbct_ig found in cbct,def.h
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

# constructor for MIRT_cbct_ig
function cbct(ig::MIRT_image_geom)
	return MIRT_cbct_ig(
		Cint.(ig.nx), Cint.(ig.ny), Cint.(ig.nz), Cfloat.(ig.dx),
		Cfloat.(ig.dy), Cfloat.(ig.dz), Cfloat.(ig.offset_x),
		Cfloat.(ig.offset_y), Cfloat.(ig.offset_z),
		pointer(ig.mask), pointer(ig.iy_start),pointer(ig.iy_end)
		)
end


# see https://docs.julialang.org/en/stable/manual/constructors/
"""
`ig = image_geom(...)`

Constructor for `MIRT_image_geom`

option:
* `nx::Integer			= 128`
* `ny::Integer			= nx`
* `dx::Real				= ?`
* `dy::Real				= -dx`
* `offset_x::Real		= 0`
* `offset_y::Real		= 0`
* `fov::Real			= ?` (if specific, then `nx*dx=ny*dy`)
* `nz::Integer			= 0`
* `dz::Real				= ?`
* `zfov::Real			= ?`
* `offset_z::Real		= 0`
* `offsets::String		= ""` (e.g., "dsp")
* `mask_type::String	= ""` (e.g., "all-but-edge")
"""
function image_geom(;
		nx::Integer			= 128,
		ny::Integer			= nx,
		dx::Real			= NaN,
		dy::Real			= NaN,
		offset_x::Real		= 0,
		offset_y::Real		= 0,
		fov::Real			= NaN,
		nz::Integer			= 0,
		dz::Real			= NaN,
		zfov::Real			= NaN,
		offset_z::Real		= 0,
		offsets::String		= "",
		mask_type::String	= "",
	) where T <: Real

	# now for the constructor
	# start by setting variables if optional arguments not used

	is3 = !isnan(dz) || !isnan(zfov) || (nz > 0)

	# offsets
	if offsets == "dsp"
		# now check for errors
		if offset_x == 0 || offset_y == 0
			throw("offsets usage incorrect")
		end
		offset_x = 0.5
		offset_y = 0.5
		if is3
			if offset_z == 0
				throw("offsets usage incorrect")
			end
		end
	end

	# transverse plane (2d) distances
	if true
		if isnan(fov)
			if isnan(dx)
				throw("dx or fov required")
			end
			fov = nx * dx
		else
			if !isnan(dx)
				throw("dx and fov?")
			end
			dx = fov / nx
			dy = fov / ny
		end

		dy = isnan(dy) ? -dx : dy # default dy
	end

	# 3D geometry
	if is3
		if isnan(zfov)
			if isnan(dz)
			#	dz = dx # nah, too risky; voxels rarely cubic
				throw("dz or zfov required")
			end
			zfov = nz * dz
		else
			if !isnan(dz)
				throw("dz and zfov?")
			end
			dz = zfov / nz
		end
	end

	# mask
	flag_mask_all_but_edge = false
	if mask_type == "all-but-edge"
		flag_mask_all_but_edge = true
		mask = trues(nx,ny,nz)
		mask[  1,   :, :] .= false
		mask[end,   :, :] .= false
		mask[  :,   1, :] .= false
		mask[  :, end, :] .= false
	end

	if mask_type == ""
		mask = is3 ? trues(nx, ny, nz) : trues(nx, ny)
	elseif size(mask,1) != nx || size(mask,2) != ny || (is3 && size(mask,3) != nz)
		throw("bad input mask size, nx=$nx ny=$ny")
	end

  # start is 0, and end is ny for a single thread
  iy_start = [0]
  iy_end   = [0]

  # image_geom_circle()
  # default is a circle that just inscribes the square
  # but keeping a 1 pixel border due to ASPIRE regularization restriction
  _circ = (;
	rx::Real = min(abs((nx/2-1)*dx), abs((ny/2-1)*dy)),
	ry::Real=rx, cx::Real=0, cy::Real=0) ->
	begin
		throw("todo: replace the 'ellipse_im' function")
		circ = ones(nx,ny)
		if is3
			circ = repmat(circ, 1, nz)
		end
		return circ
	end

  # image_geom_over()
  _over = (over::Real=1) -> MIRT_image_geom(
		nx*over, ny*over, dx/over, dy/over,
		offset_x*over, offset_y*over, fov,
		nz*over, dz/over, zfov, offset_z*over,
		offsets, mask_type, mask, iy_start, iy_end, is3,
 #todo
		_circ, _expand_nz, _over,
		_plot)

  # image_geom_plot()
  _plot = () ->
	begin
		# todo: figure out plotting with ct geom
	end

	# image_geom_expand_nz
	_expand_nz = (nz_pad::Integer=0) ->
	begin
		if !is3
			throw("expand_nz only valid for 3D")
		end
		out_nz = nz + 2*nz_pad
		out_zfov = out_nz / nz * zfov
		out_mask = cat(dims=3, repeat(mask[:,:,1], 1, 1, nz_pad), mask,
		repeat(mask[:,:,end], 1, 1, nz_pad))
		return out_nz, out_zfov, out_mask, [nx ny out_nz]
	end

	# and finally return the new object
	return MIRT_image_geom(
		# set the variables
		nx, ny, dx, dy, offset_x, offset_y, fov, nz, dz, zfov, offset_z,
		offsets, mask_type, mask, iy_start, iy_end, is3,
		# set the functions, todo
		_circ, _expand_nz, _over,
		_plot, #_wx, _wy, _wz, _u, _v, _ug, _vg, _zg, _z,
		)
end

# define a helper function needed to downsample image_geom
# image_geom_down_round()
function _down_round(val::Real=1, dd::Real=1, down::Real=1)
	out = val / down
	if out != round(out)
		out = 2 * round(out / 2)
	end
	dd = dd * down
	return Int(out), dd
end

function _downsample1(x::Array{T}, down::Integer) where {T<:Number}
	dim = size(x)
	x = reshape(x, dim[1], prod(dim[2:end]))
	m1 = floor(Int, dim[1] / down)
	if m1*down < dim[1]
		x = x[1:(m1*down), :]
	end
	y = reshape(x, down, :)
	y = mean(y,1)
	return reshape(y, (m1, dim[2], dim[3]))
end

function _downsample3(x::Array{T}, m::Array{Int}) where {T<:Number}
	ret = _downsample1(x, m[1])
	ret = _downsample1(permutedims(ret, [2, 1, 3]), m[2])
	ret = _downsample1(permutedims(ret, [3, 2, 1]), m[3])
	return permutedims(ret, [2, 1, 3])
end

# define the image_geom_downsample function as its own functions
function downsample(ig::MIRT_image_geom; down::Integer=1)

	# check if its a scalar by seeing if the collection of size is empty
	if isempty(size(down))
		if ig.is3
			downv = [ down down down]
		else
			downv = [ down down]
		end
	else
		downv = down
	end

	# now call the down round function
	down_nx, down_dx = _down_round(ig.nx, ig.dx, downv[1])
	down_ny, down_dy = _down_round(ig.ny, ig.dy, downv[2])
	# adjust to "pixel" units
	down_offset_x = ig.offset_x / downv[1]
	down_offset_y = ig.offset_y / downv[2]
	# now the 3 d case
	if ig.is3
		down_nz, down_dz = _down_round(ig.nz, ig.dz, downv[3])
		down_offset_z = ig.offset_z / downv[3]
	else
		down_nz = 0; down_dz = 0; down_offset_z = 0
	end
	# now to correctly downsample the mask
	# size is a tuple, all in Julia expects an array so we must convert
	mdim_vec = [size(ig.mask)[1] size(ig.mask)[2] size(ig.mask)[3]]
	if ig.is3
		if all([down_nx down_ny down_nz] .* downv == mdim_vec)
			down_mask = _downsample3(ig.mask, downv)
			#down_mask = trues(ig.mask)
		else
			throw("bug: bad mask size. need to address mask downsampling")
		end
	else
		if down_nx * downv[1] == size(ig.mask,1) && down_ny * downv[2] == size(ig.mask,2)
		# mask = downsample2(mask, downv) > 0
		down_mask = trues(ig.mask)
		elseif down_nx != size(ig.mask,1) || down_ny != size(ig.mask,2)
			throw("bug: bad mask size. need to address mask downsampling")
		end
	end

	return image_geom(
			# new downsampled version to be constructed
			nx			= down_nx,
			ny			= down_ny,
			dx			= down_dx,
			dy			= down_dy,
			offset_x	= down_offset_x,
			offset_y	= down_offset_y,
			nz			= down_nz,
			dz			= down_dz,
			offset_z	= down_offset_z,
			offsets		= ig.offsets,
			mask_type	= ig.mask_type,
		)
end


function Base.display(ig::MIRT_image_geom)
	ir_dump(ig)
end


# Extended properties

fun0 = Dict([
	(:help, ig -> print(image_geom_help())),

	(:dim, ig -> ig.is3 ? [ig.nx ig.ny ig.nz] : [ig.nx ig.ny]),
	(:fovs, ig -> ig.is3 ?
		[abs(ig.dx)*ig.nx, abs(ig.dy)*ig.ny, abs(ig.dz)*ig.nz] :
			[abs(ig.dx)*ig.nx, abs(ig.dy)*ig.ny]),

	(:zeros, ig -> zeros(Float32, ig.dim...)), # (nx, ny, nz) : zeros(nx, ny)
	(:ones, ig -> ones(Float32, ig.dim...)), # (nx, ny, nz) : zeros(nx, ny)

# _ones = () -> is3 ? ones(nx, ny, nz) : ones(nx, ny)

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

	(:mask_or, ig -> ig.is3 ? dropdims(sum(ig.mask, dims=3) .> 0, dims=3) :
							ig.mask),
	(:np, ig -> sum(ig.mask)),

	(:mask_outline, ig ->
		begin
			mask2 = ig.mask_or
			tmp = conv2(Float32.(mask2), ones(Float32,3,3)) # todo: imfilter
			tmp = tmp[2:end-1,2:end-1] # 'same'
			return (tmp .> 1) .& (.! mask2)
		end),

	# functions

	(:embed, ig -> (x::AbstractArray{<:Number} -> embed(x, ig.mask))),

	(:maskit, ig -> (x::AbstractArray{<:Number} -> maskit(x, ig.mask))),

	(:shape, ig -> (x::AbstractArray{<:Number} -> reshape(x, ig.dim...))),

	(:unitv, ig -> ((;kwargs...) -> image_geom_add_unitv(ig.zeros; kwargs...))),

	])


# Tricky overloading here!

Base.getproperty(ig::MIRT_image_geom, s::Symbol) =
		haskey(fun0, s) ? fun0[s](ig) :
		getfield(ig, s)

Base.propertynames(ig::MIRT_image_geom) =
	(fieldnames(typeof(ig))..., keys(fun0)...)


function image_geom_test2()
	ig = image_geom(nx=16, dx=2)
	image_geom_test2(ig)
	true
end

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
	@test ig.embed(ig.ones[ig.mask]) == Float32.(ig.mask)
	@test ig.maskit(ig.ones) == Float32.(ones(ig.np))
	ig.mask_outline
	@test ig.shape(ig.ones[:]) == ig.ones
	ig.unitv(j=4) 
	ig.unitv(i=ones(Int, length(ig.dim)))
	ig.unitv(c=zeros(Int, length(ig.dim)))
	ig.unitv() 
	ig.ones
	ig.zeros
#	ig.circ()
	ig.over(2)
#	ig.plot()
	ig.u
	ig.v
	ig.ug
	ig.vg
	ig.fg
	true
end

function image_geom_test3(ig::MIRT_image_geom)
	image_geom_test2(ig)
	ig.wz
	ig.zg
	ig.expand_nz(2)
	ig.mask_or
	true
end

function image_geom_test3()
	ig = image_geom(nx=16, nz=3, dx=2, dz=4)
	image_geom_test3(ig)
	true
end


"""
`image_geom(:test)`
"""
function image_geom(test::Symbol)
	if test == :help
		print(image_geom_help())
		return true
	end
	test != :test && throw(ArgumentError("test $test"))
	@test image_geom_test2()
	@test image_geom_test3()
	true
end

