#=
image_geom.jl
Methods related to an image geometry for image reconstruction
2017-10-02 Samuel Rohrer, University of Michigan
2019-03-05 Jeff Fessler, Julia 1.1 + tests
2019-06-23 Jeff Fessler, overhaul
2020-07-21 Jeff Fessler, D-dimensional
=#

export ImageGeom, image_geom, cbct
export image_geom_circle
export image_geom_ellipse

#using MIRT: downsample2, downsample3
using FillArrays: Trues
using ImageTransformations: imresize
#import Unitful # Length


"""
    ImageGeom{D,S}

- `dims::Dims{D}` image dimensions
- `deltas::S where S <: NTuple{D,Number}` pixel sizes,
  where each `Number` is usually `Real` or `Unitful.Length`
- `offsets::NTuple{D,Float32}` unitless
- `mask::AbstractArray{Bool,D}` logical mask

Image geometry "struct" with essential parameters
"""
struct ImageGeom{D, S <: NTuple{D,Number}}
    dims::Dims{D} # image dimensions
    deltas::S # pixel sizes
    offsets::NTuple{D,Float32} # unitless
    mask::AbstractArray{Bool,D} # logical mask

    function ImageGeom{D,S}(
        dims::Dims{D},
        deltas::S,
        offsets::NTuple{D,Real},
        mask::AbstractArray{Bool,D},
    ) where {D, S <: NTuple{D,Number}} # <: Union{Real,Unitful.Length}
        any(x -> x <= zero(x), dims) && throw("dims must be positive")
        any(iszero, deltas) && throw("deltas must be nonzero")
        size(mask) == dims ||
            throw(DimensionMismatch("mask size $(size(mask)) vs dims $dims"))
        new{D,S}(dims, deltas, Float32.(offsets), mask)
    end
end

"""
    ig = ImageGeom(dims, deltas, offsets, [, mask])
Convenient constructor for `ImageGeom`.
The `deltas` elements should each be `Real` or a `Unitful.Length`.
Default `mask` is `FillArrays.Trues(dims)` which is akin to `trues(dims)`.
"""
function ImageGeom(
    dims::Dims{D},
    deltas::S,
    offsets::NTuple{D,Real},
    mask::AbstractArray{Bool,D},
) where {D, S <: NTuple{D,Number}}
   ImageGeom{D,S}(dims, deltas, offsets, mask)
end

function ImageGeom(
    dims::Dims{D},
    deltas::S,
    offsets::NTuple{D,Real},
) where {D, S <: NTuple{D,Number}}
   ImageGeom(dims, deltas, offsets, Trues(dims))
end


"""
image_geom_help( ; io)
"""
function image_geom_help( ; io::IO = isinteractive() ? stdout : devnull)
	print(io, "propertynames:\n")
	print(io, propertynames(image_geom(nx=1, dx=1)))

	print(io,
	"
	Derived values for 2D (and 3D) images:

	is3	is it 3d (nz > 0?)
	fovs (|dx|*nx, |dy|*ny, ...)
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
	unitv(...)	j=j | i=(ix,iy) | c=(cx cy)
				j: single index from 1 to length(z)
				i: (ix,iy[,iz]) index from 1 to nx,ny
				c: (cx,cy[,cz]) index from +/- n/2 center at floor(n/2)+1
	circ(rx=,ry=,cx=,cy=)	circle of given radius and center (cylinder in 3D)
	plot(jim)	plot the image geometry using the `MIRTjim.jim` function

	Methods that return a new `ImageGeom:`

	down(down::Int)		down-sample geometry by given factor
	over(over::Int)		over-sample geometry by given factor
	expand_nz(nz_pad)	expand image geometry in z by nz_pad on both ends
	\n")
end


"""
    MIRT_cbct_ig

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
    cbct(ig::ImageGeom{3,<:Real}; nthread::Int=1)
Constructor for `MIRT_cbct_ig` (does not support units currently)
"""
function cbct(ig::ImageGeom{3,S} ; nthread::Int=1) where {S <: NTuple{3,Real}}
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


"""
    ig = image_geom(...)

Constructor for `ImageGeom`, where `dx,dy,dz` and `fov` and `fovz` may have units

# Arguments
- `nx::Int = 128`
- `ny::Int = nx`
- `dx::Number = ?` (must specify one of `dx` or `fov`)
- `dy::Number = -dx`
- `offset_x::Real = 0` (unitless)
- `offset_y::Real = 0` (unitless)
- `fov::Number = ?` (if specified, then `nx*dx=ny*dy`)
- `nz::Int = 0`
- `dz::Number = ?` (need one of `dz` or `zfov` if `nz > 0`)
- `zfov::Number = ?` (if specified, then `nz*dz`)
- `offset_z::Real = 0` (unitless)
- `offsets::Symbol = :none` or :dsp
- `mask::Union{Symbol,AbstractArray{Bool}} = :all` | `:circ` | `:all_but_edge_xy`
"""
function image_geom( ;
	nx::Int = 128,
	ny::Int = nx,
	dx::Number = NaN,
	dy::Number = NaN,
	offset_x::Real = 0,
	offset_y::Real = 0,
	fov::Number = NaN,
	nz::Int = 0,
	dz::Number = NaN,
	zfov::Number = NaN,
	offset_z::Real = 0,
	offsets::Symbol = :none,
	mask::Union{Symbol,AbstractArray{Bool}}	= :all,
)

	# handle optional arguments

	is3 = !isnan(dz) || !isnan(zfov) || (nz > 0)

	# offsets
	if offsets === :dsp
		# check that no nonzero offset was specified
		(offset_x != 0 || offset_y != 0) && throw("offsets usage incorrect")
		offset_x = 0.5
		offset_y = 0.5
		if is3
			offset_z != 0 && throw("offset_z usage incorrect")
			offset_z = 0.5
		end
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
	if mask === :all
		mask = is3 ? Trues(nx, ny, nz) : Trues(nx, ny)
	elseif mask === :circ
		mask = image_geom_circle(nx, ny, dx, dy)
		if nz > 0
			mask = repeat(mask, 1, 1, nz)
		end
	elseif mask === :all_but_edge_xy
		mask = trues(nx,ny,max(nz,1))
		mask[  1,   :, :] .= false
		mask[end,   :, :] .= false
		mask[  :,   1, :] .= false
		mask[  :, end, :] .= false
		if !is3
			mask = mask[:,:,1]
		end
	elseif isa(mask, Symbol)
		throw("mask symbol $mask")
	elseif size(mask,1) != nx || size(mask,2) != ny || (is3 && size(mask,3) != nz)
		throw("mask size $(size(mask)), nx=$nx ny=$ny nz=$nz")
	end

	# return the object (type unstable, but used rarely so ok)
	return is3 ?
		ImageGeom((nx, ny, nz), (dx, dy, dz),
			(offset_x, offset_y, offset_z), mask) :
		ImageGeom((nx, ny), (dx, dy),
			(offset_x, offset_y), mask)
end



"""
    (n,d) = _down_round()
helper function needed to downsample `image_geom`
"""
function _down_round(
	val::NTuple{D,Real},
	dd::NTuple{D,Number},
	down::NTuple{D,Int},
) where {D}
	# for non-divisors make dim a multiple of 2
	fun = r -> r == round(r) ? Int(r) : 2 * round(Int, r/2)
	return @. fun(val / down), dd .* down
end


"""
    ig_down = ig_downsample(ig, down::Tuple{Int})
down sample an image geometry by the factor `down`
cf `image_geom_downsample`
"""
function ig_downsample(ig::ImageGeom{D,S}, down::NTuple{D,Int}) where {D,S}

	# call the down round function
	down_dim, deltas = _down_round(ig.dims, ig.deltas, down)
	# adjust offsets to new "pixel" units
	down_offsets = ig.offsets ./ down

	# carefully down-sample the mask
	mdim = size(ig.mask)
	if ig.is3
		if down_dim .* down == mdim
			down_mask = downsample3(ig.mask, down) .> 0
		else
			throw("bug: bad mask size. need to address mask downsampling")
		end
	else
		if down_dim .* down == mdim
			down_mask = downsample2(ig.mask, down) .> 0
		else
			down_mask = imresize(ig.mask, down_dim...) .> 0
		end
	end

	return ImageGeom{D,S}(down_dim, deltas, down_offsets, down_mask)
end

ig_downsample(ig::ImageGeom{D}, down::Int) where {D} =
	ig_downsample(ig, ntuple(i->down, D))


"""
ig_new = image_geom_expand_nz(ig::ImageGeom{3}, nz_pad::Int)
pad both ends
"""
function image_geom_expand_nz(ig::ImageGeom{3,S}, nz_pad::Int) where S
	out_nz = ig.nz + 2*nz_pad
	out_mask = cat(dims=3, repeat(ig.mask[:,:,1], 1, 1, nz_pad),
		ig.mask, repeat(ig.mask[:,:,end], 1, 1, nz_pad),
	)
	return ImageGeom{3,S}((ig.dims[1], ig.dims[2], out_nz),
		ig.deltas, ig.offsets, out_mask,
	)
end


"""
ig_over = image_geom_over(ig::ImageGeom, over::Int)
over-sample an image geometry by the factor `over`
"""
function image_geom_over(ig::ImageGeom{D}, over::Int) where {D}
	if ig.mask == Trues(ig.dims) || all(ig.mask)
		mask_over = Trues(ig.dims .* over)
	else
		mask_over = imresize(ig.mask, ig.dims .* over) .> 0
	end
	return ImageGeom(ig.dims .* over, ig.deltas ./ over,
		ig.offsets .* over, mask_over,
	)
end


# ellipse that just inscribes the rectangle
# but keeping a 1 pixel border due to ASPIRE regularization restriction
function image_geom_ellipse(
	nx::Int, ny::Int, dx::Number, dy::Number ;
	rx::Number = min(abs((nx/2-1)*dx), abs((ny/2-1)*dy)),
	ry::Number = min(abs((nx/2-1)*dx), abs((ny/2-1)*dy)),
	cx::Number = 0, cy::Number = 0, over::Int=2,
)
	ig = image_geom(nx=nx, ny=ny, dx=dx, dy=dy)
	circ = ellipse_im(ig, [cx cy rx ry 0 1], oversample=over) .> 0
	return circ
end


# default is a circle that just inscribes the square
# but keeping a 1 pixel border due to ASPIRE regularization restriction
function image_geom_circle(
	nx::Int, ny::Int, dx::Number, dy::Number ;
	rx::Number = min(abs((nx/2-1)*dx), abs((ny/2-1)*dy)),
	ry::Number = rx, cx::Number = 0, cy::Number = 0, nz::Int=0, over::Int=2,
)
	ig = image_geom(nx=nx, ny=ny, dx=dx, dy=dy)
	circ = ellipse_im(ig, [cx cy rx ry 0 1], oversample=over) .> 0
	if nz > 0
		circ = repeat(circ, 1, 1, nz)
	end
	return circ
end


"""
out = image_geom_add_unitv(z::AbstractArray ; j=?, i=?, c=?)

add a unit vector to an initial array `z` (typically of zeros)

# options (use at one of these):
- `j` single index from 1 to length(z)
- `i` (ix, iy [,iz]) index from 1 to nx,ny
- `c` (cx, cy [,cz]) index from +/- n/2 center at floor(n/2)+1

default with no arguments gives unit vector at center `c=(0, 0 [,0])`
"""
function image_geom_add_unitv(
	z::AbstractArray{T, D} ; # starts with zeros()
	j::Int = 0,
	i::NTuple{D,Int} = ntuple(i -> 0, D),
	c::NTuple{D,Int} = ntuple(i -> 0, D),
) where {T <: Number, D}

	out = copy(z)

	(j < 0 || j > length(z)) && throw("bad j $j")
	if 1 <= j <= length(z)
		any(!=(0), i) && throw("i $i")
		any(!=(0), c) && throw("c $c")
		out[j] += one(T)
	elseif all(1 .<= i .<= size(z))
		any(!=(0), c) && throw("c $c")
		out[i...] += one(T)
	else
		any(!=(0), i) && throw("i $i")
		tmp = c .+ (size(out) .÷ 2) .+ 1
		out[tmp...] += one(T)
	end

	return out
end


"""
    image_geom_plot(ig, how ; kwargs...)
The `how` argument should be `MIRTjim.jim` to be useful.
"""
image_geom_plot(ig::ImageGeom{2}, how::Function ; kwargs...) =
	how(ig.x, ig.y, ig.mask, "(nx,ny)=$(ig.nx),$(ig.ny)" ; kwargs...)
image_geom_plot(ig::ImageGeom{3}, how::Function ; kwargs...) =
	how(ig.x, ig.y, ig.mask_or,
		"(dx,dy,dz)=$(ig.dx),$(ig.dy),$(ig.dz)" ; kwargs...)


"""
    show(io::IO, ig::ImageGeom)
    show(io::IO, ::MIME"text/plain", ig::ImageGeom)
"""
Base.show(io::IO, ig::ImageGeom) =
	print(io, "ImageGeom: $(ig.dim)")

function Base.show(io::IO, ::MIME"text/plain", ig::ImageGeom)
	ir_dump(io, ig)
end

_zero(ig::ImageGeom{D,S}) where {D, S <: NTuple{D,Real}} = zero(Float32)
_zero(ig::ImageGeom{D,S}) where {D, S <: NTuple{D,Any}} = zero(Int32) # alert!
_zero(ig::ImageGeom{D,S}) where {D, S <: NTuple{D,T}} where {T <: Number} = zero(T)
_zero(ig::ImageGeom{D,S}) where {D, S <: NTuple{D,T}} where {T <: Real} = zero(T)

# Extended properties

image_geom_fun0 = Dict([
	(:help, ig -> image_geom_help()),

	(:ndim, ig -> length(ig.dims)),
	(:nx, ig -> ig.dims[1]),
	(:ny, ig -> ig.ndim >= 2 ? ig.dims[2] : 0),
	(:nz, ig -> ig.ndim >= 3 ? ig.dims[3] : 0),
	(:is3, ig -> ig.nz > 0),
	(:dim, ig -> ig.dims),
	(:fovs, ig -> abs.(ig.deltas) .* ig.dims),

	(:zeros, ig -> zeros(Float32, ig.dim)),
	(:ones, ig -> ones(Float32, ig.dim)),

	(:dx, ig -> ig.deltas[1]),
    (:dy, ig -> ig.ndim >= 2 ? ig.deltas[2] : _zero(ig)),
    (:dz, ig -> ig.ndim >= 3 ? ig.deltas[3] : _zero(ig)),

	(:offset_x, ig -> ig.offsets[1]),
	(:offset_y, ig -> ig.ndim >= 2 ? ig.offsets[2] : Float32(0)),
	(:offset_z, ig -> ig.ndim >= 3 ? ig.offsets[3] : Float32(0)),

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
						zeros(Float32, ig.nx, ig.ny)),

	# DFT frequency sample grid
	(:u, ig -> (-ig.nx/2:(ig.nx/2-1)) / (ig.nx*ig.dx)),
	(:v, ig -> (-ig.ny/2:(ig.ny/2-1)) / (ig.ny*ig.dy)),
	(:w, ig -> (-ig.nz/2:(ig.nz/2-1)) / (ig.nz*ig.dz)),

	(:ug, ig -> ig.is3 ? repeat(ig.u, 1, ig.ny, ig.nz) :
						repeat(ig.u, 1, ig.ny)),
	(:vg, ig -> ig.is3 ? repeat(ig.v', ig.nx, 1, ig.nz) :
						repeat(ig.v', ig.nx, 1)),
	(:wg, ig -> ig.is3 ? repeat(reshape(ig.w, 1, 1, ig.nz), ig.nx, ig.ny, 1) :
						zeros(Float32, ig.nx, ig.ny)),
	(:fg, ig -> ig.is3 ? (ig.ug, ig.vg, ig.wg) : (ig.ug, ig.vg)),

	(:np, ig -> sum(ig.mask)),
	(:mask_or, ig -> mask_or(ig.mask)),
	(:mask_outline, ig -> mask_outline(ig.mask)),

	# simple functions

	(:plot, ig -> ((how::Function ; kwargs...) ->
		image_geom_plot(ig, how ; kwargs...))),
	(:embed, ig -> (x::AbstractArray{<:Number} -> embed(x, ig.mask))),
	(:maskit, ig -> (x::AbstractArray{<:Number} -> maskit(x, ig.mask))),
	(:shape, ig -> (x::AbstractArray{<:Number} -> reshape(x, ig.dim))),
	(:unitv, ig -> (( ; kwargs...) -> image_geom_add_unitv(ig.zeros ; kwargs...))),
	(:circ, ig -> (( ; kwargs...) ->
		image_geom_circle(ig.nx, ig.ny, ig.dx, ig.dy, nz=ig.nz ; kwargs...))),

	# functions that return new geometry:

	(:down, ig -> (down::Int -> ig_downsample(ig, down))),
	(:over, ig -> (over::Int -> image_geom_over(ig, over))),
	(:expand_nz, ig -> (nz_pad::Int -> image_geom_expand_nz(ig, nz_pad))),

	])


# Tricky overloading here!

Base.getproperty(ig::ImageGeom, s::Symbol) =
		haskey(image_geom_fun0, s) ? image_geom_fun0[s](ig) :
		getfield(ig, s)

Base.propertynames(ig::ImageGeom) =
	(fieldnames(typeof(ig))..., keys(image_geom_fun0)...)
