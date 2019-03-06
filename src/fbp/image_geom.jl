# image_geom.jl
# Methods related to an image geometry for image reconstruction
# 2017-10-02 Samuel Rohrer, University of Michigan
# 2019-03-05 Jeff Fessler, Julia 1.1 + tests

using Printf
using Test
using DSP: conv2
using SparseArrays: issparse#, find, ind,

"""
`array = embed(v, mask)`

embed vector v of length sum(mask) into elements of an array where mask is true
"""
function embed(v::AbstractArray{T,1} where T <: Number,
		mask::AbstractArray{Bool,N} where N)
	array = zeros(eltype(v), size(mask))
	array[mask] .= v
	return array
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
		i::AbstractVector{T1} = zeros(Int, ndims(z)),
		c::AbstractVector{T2} = zeros(Int, ndims(z))
		) where {T1 <: Integer, T2 <: Integer}
	out = copy(z)

	if j > 0 && all(i .== 0)
		out[j] += 1
	elseif j == 0 && any(i .> 0)
		out[i...] += 1
		#	elseif is3 & ix > 0 && iy > 0 && iz > 0
		#		out[ix,iy,iz] .= 1
		#	elseif !is3 & ix > 0 && iy > 0
		#		out[ix,iy] .= 1
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
	@assert test == :test
	image_geom_add_unitv(zeros(3,4), j=2)
	image_geom_add_unitv(zeros(3,4), i=[2,3])
	image_geom_add_unitv(zeros(3,4))
	image_geom_add_unitv(zeros(3,4), c=[-1,-2])
	true
end


# Image geometry "struct" with some parameters and methods included
struct MIRT_image_geom
  # options for 2D image geometry
  nx::Integer           # image dimension one
  ny::Integer           # image dimension (default: nx) two
  dx::Real              # pixel size ('dx' or 'fov' required)
  dy::Real              # pixel size (default: -dx). (value | 'dx' | '-dx' | '2*dx')
  offset_x::Real        # unitless (default: 0)
  offset_y::Real        # unitless (default: 0)
  fov::Real             # nx*dx = ny*dy (if specified)

  # options for 3D image geometry
  nz::Integer           # image dimension three
  dz::Real              # pixel size (default: dx)
  zfov::Real            # dz*nz
  offset_z::Real        # unitless (default: 0)

  # general image orientation
  offsets::String       # this is a string, only 'dsp' has meaning now
  mask_type::String     # a string stored to make downsample easier
  mask::Array{Bool}     # logical mask
  iy_start::Array{Int}  # size [nthread]
  iy_end::Array{Int}    # size [nthread]
  is3::Bool             # storing this to make life easier

  # functions provided by the constructor
  # methods for 2D images
  dim::Function              # [nx ny [nz]]
  x::Function                # (subs) 1D x coordinates of each pixel
  y::Function                # (subs) 1D y coordiantes of each pixel
  fg::Function               # {subs} cell of 2D or 3D frequency coordinates
  xg::Function               # x coordinates of each pixel as a grid (2D or 3D)
  yg::Function               # y coordinates of each pixel as a grid (2D or 3D)
  fovs::Function             # [|dx|*nx |dy|*ny ...]
  np::Function               # sum(mask(:)) = # pixels to be estimated
  embed::Function            # turn short column(s) into array(s)
  maskit::Function           # opposite of embed
  mask_outline::Function     # binary image showing 2D outline of mask
  shape::Function            # reshape long column x to [nx ny [nz]] array
  unitv::Function            # (ix,iy) | (jj) | ('c', [cx cy])
  ones::Function             # | ('nz', nz) | (jj)
  zeros::Function            # | ('nz', nz) | (jj)
  circ::Function             # (rx,ry,cx,cy) circle of radius rx, ry (cylinder in 3D)
  expand_nz::Function        # '(nx_pad)'
  over::Function             #
  plot::Function             # plot the image geometry
  wx::Function               # '(), (nx - 1/2) * dx + offset_x'
  wy::Function               # '(), (ny - 1/2) * dy + offset_y'
  wz::Function               # '(), (nz - 1/2) * dz + offset_z'
  u::Function                # 1D frequency domain coordinates of each pixel
  v::Function                # 1D frequency domain coordinates of each pixel
  ug::Function               # 2D or 3D grid of frequency domain coordinates
  vg::Function               # 2D or 3D grid of frequency domain coordinates

  # methods for 3D images
  zg::Function              # z coordinates of each pixel as a grid (3D only)
  z::Function               # z coordinates of each pixel
  mask_or::Function         # [nx ny] logical "or" of mask
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
Constructor for MIRT_image_geom
ig = image_geom(...)

option:
	nx::Integer              = 128
	ny::Integer              = nx
	dx::Real                 = ?
	dy::Real                 = -dx
	offset_x::Real           = 0
	offset_y::Real           = 0
	fov::Real                = ? (if specific, then nx*dx=ny*dy)
   	nz::Integer              = 0
   	dz::Real                 = ?
   	zfov::Real               = ?
   	offset_z::Real           = 0
   	offsets::String          = "" (e.g., "dsp")
   	mask_type::String        = "" (e.g., "all-but-edge")
"""
function image_geom(;
		nx::Integer              = 128,
		ny::Integer              = nx,
		dx::Real                 = NaN,
		dy::Real                 = NaN,
		offset_x::Real           = 0,
		offset_y::Real           = 0,
		fov::Real				= NaN,
		nz::Integer              = 0,
		dz::Real                 = NaN,
		zfov::Real               = NaN,
		offset_z::Real           = 0,
		offsets::String          = "",
		mask_type::String        = "",
	) where T <: Real

  # now for the constructor
  # start by setting variables if optional arguments not used

  is3 = !isnan(dz) || !isnan(zfov) || (nz > 0)

  # offsets
  if offsets == "dsp"
    # now check for errors
    if offset_x == 0 || offset_y == 0
      error("offsets usage incorrect")
    end
    offset_x = 0.5
    offset_y = 0.5
    if is3
      if offset_z == 0
        error("offsets usage incorrect")
      end
    end
  end

	# transverse plane (2d) distances
	if true
		if isnan(fov)
			if isnan(dx)
				error("dx or fov required")
			end
			fov = nx * dx
		else
			if !isnan(dx)
				error("dx and fov?")
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
				error("dz or zfov required")
			end
			zfov = nz * dz
		else
        	if !isnan(dz)
				error("dz and zfov?")
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
    error(@sprintf("bad input mask size, nx=%d ny=%d", nx, ny))
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
		error("todo: replace the 'ellipse_im' function")
		circ = ones(nx,ny)
		if is3
			circ = repmat(circ, 1, nz)
		end
		return circ
	end

  # _dim
  _dim = () -> return is3 ? [nx ny nz] : [nx ny]

  # image_geom_mask_or()
  _mask_or = () -> reshape(sum(mask, dims=3) .> 0, nx,ny) # squeeze??

  # image_geom_np()
  _np = () -> sum(mask[:])

  # image_geom_over()
  _over = (over::Real=1) -> MIRT_image_geom(
        nx*over, ny*over, dx/over, dy/over,
		offset_x*over, offset_y*over, fov,
		nz*over, dz/over, zfov, offset_z*over,
        offsets, mask_type, mask, iy_start, iy_end, is3,
        _dim, _x, _y, _fg, _xg, _yg, _fovs, _np, _embed, _maskit, _mask_outline,
        _shape, _unitv, _ones, _zeros, _circ, _expand_nz, _over,
        _plot, _wx, _wy, _wz, _u, _v, _ug, _vg, _zg, _z, _mask_or,)

  # image_geom_embed_sparse
  # called by image_geom_embed
  function _embed_sparse(x::Array{T}) where {T <: Number}
      i, j, a = find(x)
      ind = find(mask)
      j = ind(j)
      return x = sparse(i, j, a, size(x,1), prod(_dim()))
  end

  # image_geom_embed
  function _embed(x::Array{T}) where {T <: Number}
	return embed(issparse(x) ? _embed_sparse(x) : x, mask)
  end

  # image_geom_maskit()
  function _maskit(x::Array{T}) where {T <: Number}
      dim = size(x)
      x = reshape(x, prod(_dim()))
      x = x[reshape(mask, prod(_dim()))]
      if length(dim) > length(_dim())
        x = reshape(x, dim[(1+length(_dim())):end])
      end
      return x
    end

  # image_geom_mask_outline()
  _mask_outline = () ->
	begin
		mask2 = is3 ? _mask_or() : mask
		tmp = conv2(Float32.(mask2), ones(Float32,3,3))
		tmp = tmp[2:end-1,2:end-1] # 'same'
		return (tmp .> 1) .& (.! mask2)
    end

  # image_geom_plot()
  _plot = () ->
    begin
      # todo: figure out plotting with ct geom
    end

  # image_geom_shape()
  function _shape(x::Array{T}) where {T <: Number}
      return is3 ? reshape(x, nx, ny, nz) : reshape(x, nx, ny)
  end

  # image_geom_wx()
  _wx = () -> (nx - 1)/2 + offset_x

  # image_geom_wy()
  _wy = () -> (ny - 1)/2 + offset_y

  # image_geom_wz()
  _wz = () -> (nz - 1)/2 + offset_z

  # image_geom_x()
  #_x = () -> (collect(0:nx-1) .- _wx()) * dx # todo
  _x = () -> ((0:nx-1) .- _wx()) * dx

  # image_geom_y()
  #_y = () -> (collect(0:ny-1) .- _wy()) * dy
  _y = () -> ((0:ny-1) .- _wy()) * dy

  # image_geom_z()
  #_z = () -> (collect(0:nz-1) .- _wz()) * dz
  _z = () -> ((0:nz-1) .- _wz()) * dz

  # image_geom_xg()
  _xg = () -> is3 ? repeat(_x(), 1, ny, nz) : repeat(_x(), 1, ny)

  # image_geom_yg()
  _yg = () -> is3 ? repeat(_y()', nx, 1, nz) : repeat(_y()', nx, 1)

  # image_geom_zg()
  _zg = () -> is3 ? repeat(reshape(_z(), 1, 1, nz), nx, ny, 1) : zeros(nx,ny)

  # image_geom_u()
  _u = () -> (-nx/2:nx/2-1) / (nx*dx)

  # image_geom_v()
  _v = () -> (-ny/2:ny/2-1) / (ny*dy)

  # image_geom_fovs
  _fovs = () -> is3 ? [abs(dx)*nx, abs(dy)*ny, abs(dz)*nz] : [abs(dx)*nx, abs(dy)*ny]

  # DFT frequency sample grid

  # image_geom_ug
  _ug = () -> is3 ? repeat(_u(), 1, ny, nz) : repeat(_u(), 1, ny)

  # image_geom_vg
  _vg = () -> is3 ? repeat(_v()', nx, 1, nz) : repeat(_v()', nx, 1)

  _wg = () -> is3 ? repeat(reshape((-nz/2:nz/2-1)/nz/dz, 1, 1, nz), nx, ny, 1) : zeros(nx,ny)

  # image_geom_fg()
  _fg = () -> is3 ? [_ug(), _vg(), _wg()] : [_ug(), _vg()]

  # image_geom_zeros
  _zeros = () -> is3 ? zeros(nx, ny, nz) : zeros(nx, ny)

  # image_geom_ones
  _ones = () -> is3 ? ones(nx, ny, nz) : ones(nx, ny)

  # image_geom_expand_nz
  _expand_nz = (nz_pad::Integer=0) ->
    begin
      if !is3
        error("expand_nz only valid for 3D")
      end
      out_nz = nz + 2*nz_pad
      out_zfov = out_nz / nz * zfov
      out_mask = cat(dims=3, repeat(mask[:,:,1], 1, 1, nz_pad), mask,
      		repeat(mask[:,:,end], 1, 1, nz_pad))
      return out_nz, out_zfov, out_mask, [nx ny out_nz]
    end

  # image_geom_unitv
  _unitv = (;kwargs...) -> image_geom_add_unitv(_zeros(); kwargs...)

  # and finally return the new object
  return MIRT_image_geom(
        # set the variables
        nx, ny, dx, dy, offset_x, offset_y, fov, nz, dz, zfov, offset_z,
        offsets, mask_type, mask, iy_start, iy_end, is3,
        # set the functions
        _dim, _x, _y, _fg, _xg, _yg, _fovs, _np, _embed, _maskit, _mask_outline,
        _shape, _unitv, _ones, _zeros, _circ, _expand_nz, _over,
        _plot, _wx, _wy, _wz, _u, _v, _ug, _vg, _zg, _z, _mask_or,)
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
  y = reshape(y, (m1, dim[2], dim[3]))
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
      error("bug: bad mask size. need to address mask downsampling")
    end
  else
    if down_nx * downv[1] == size(ig.mask,1) && down_ny * downv[2] == size(ig.mask,2)
      # mask = downsample2(mask, downv) > 0
      down_mask = trues(ig.mask)
    elseif down_nx != size(ig.mask,1) || down_ny != size(ig.mask,2)
      error("bug: bad mask size. need to address mask downsampling")
    end
 end

  return image_geom(
           # new downsampled version to be constructed
           nx         = down_nx,
           ny         = down_ny,
           dx         = down_dx,
           dy         = down_dy,
           offset_x   = down_offset_x,
           offset_y   = down_offset_y,
           nz         = down_nz,
           dz         = down_dz,
           offset_z   = down_offset_z,
           offsets    = ig.offsets,
           mask_type  = ig.mask_type,
         )
end


function Base.display(ig::MIRT_image_geom)
	ir_dump(ig)
end


function image_geom_test2()
	ig = image_geom(nx=16, dx=2)
	image_geom_test2(ig)
	true
end

function image_geom_test2(ig::MIRT_image_geom)
	# test 2D functions provided by the constructor
	ig.dim()
	ig.x()
	ig.y()
	ig.fg()
	ig.xg()
	ig.yg()
	ig.fovs()
	ig.np()
	ig.embed(ig.ones()[ig.mask])
	ig.maskit(ig.ones())
	ig.mask_outline()
	ig.shape(ig.ones()[:])
	ig.unitv(j=4) 
	ig.unitv(i=ones(Int, length(ig.dim())))
	ig.unitv(c=zeros(Int, length(ig.dim())))
	ig.unitv() 
	ig.ones()
	ig.zeros()
#	ig.circ()
	ig.over(2)
#	ig.plot()
	ig.wx()
	ig.wy()
	ig.wz()
	ig.u()
	ig.v()
	ig.ug()
	ig.vg()
	true
end

function image_geom_test3(ig::MIRT_image_geom)
	image_geom_test2(ig)
	ig.expand_nz(2)
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
	@assert test == :test
	@test image_geom_test2()
	@test image_geom_test3()
	true
end


"""
`embed(:test)`
"""
function embed(test::Symbol)
	@assert test == :test
	mask = [false true true; true false false]
	@assert embed(1:3,mask) == [0 2 3; 1 0 0]
	true
end
