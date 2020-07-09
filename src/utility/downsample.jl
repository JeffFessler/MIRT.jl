#=
downsample.jl
Copyright 2019-03-05, Jeff Fessler, University of Michigan

see also MIRT/time/downsample.jl
=#

export downsample_dim1
export downsample1
export downsample2
export downsample3
export downsample


"""
    y = downsample_dim1(x, down ; warn::Bool)
downsample by factor `down` along first dimension by averaging

in
- `x [n1 (Nd)]`
- `down::Int` downsampling factor

option
- `warn::Bool` warn if noninteger multiple; default `isinteractive()`

out
- `y [n1/down (Nd)]`

Copyright 2019-03-05, Jeff Fessler, University of Michigan
"""
function downsample_dim1(x::AbstractArray{<:Number}, down::Int
		; warn::Bool = isinteractive())

	dim = size(x)
	dim1 = dim[1]
	x = reshape(x, dim1, :) # [n1 *Nd]
	m1 = Int(floor(dim1 / down))
	if m1 * down < dim1
		warn && @warn("truncating input size $dim1 to $(m1 * down)")
		x = x[1:(m1*down),:]
	end
	y = reshape(x, down, :)
	y = sum(y, dims=1) / down # mean
	y = reshape(y, m1, dim[2:end]...)
#	out = similar(x, (m1, dim[2:end]...))
#	vec(out) .= vec(y) # failed attempt to help @inferred for 2D input arrays
	return y
end


"""
    y = downsample1(x, down ; warn=true)
downsample 1D vector by factor `down`

in
- `x [n1]`
- `down::Int` downsampling factor

option
- `warn::Bool` warn if noninteger multiple; default `isinteractive()`

out
- `y [n1/down]`

Copyright 2019-03-05, Jeff Fessler, University of Michigan
"""
function downsample1(x::AbstractVector{<:Number}, down::Int
		; warn::Bool = isinteractive())

	dim = size(x)
	dim1 = dim[1]
	m1 = floor(Int, dim1 / down)
	if m1 * down < dim1
		warn && @warn("truncating input size $dim1 to $(m1 * down)")
		y = reshape((@view x[1:(m1*down)]), down, :)
	else
		y = reshape(x, down, :)
	end
	y = sum(y, dims=1) / down # mean
	y = reshape(y, m1, dim[2:end]...)
	return vec(y)
end


"""
    y = downsample2(x, down ; warn=true, T)

downsample by averaging by integer factors
in
- `x [nx ny]`
- `down` can be a scalar (same factor for both dimensions) or a 2-vector

option
- `warn::Bool` warn if noninteger multiple; default `isinteractive()`
- `T::DataType` specify output eltype; default `eltype(x[1] / down[1])`

out
- `y [nx/down ny/down]`
"""
function downsample2(x::AbstractMatrix{<:Number},
		down::AbstractVector{Int},
		; warn::Bool = isinteractive(),
		T::DataType = eltype(x[1] / down[1])
	)

	length(down) != 2 && throw("bad down $down")

	idim = size(x)
	odim = floor.(Int, idim ./ down)

	warn && any(odim .* down .!= idim) && @warn("truncating to $odim")

#	y = similar(x, odim) # fails!?
	y = Array{T}(undef, odim[1], odim[2])
	d1 = down[1]
	d2 = down[2]
	for i2=1:odim[2]
		for i1=1:odim[1]
			y[i1,i2] =
				sum(@view x[(i1-1)*d1 .+ (1:d1), (i2-1)*d2 .+ (1:d2)]) / d1 / d2
		end
	end

#=
	fun = (x, d) -> downsample_dim1(x, d, warn=warn)

#=
#	this old way returned an Adjoint type:
	y = fun(x, down[1])
	y = fun(y', down[end])'
=#

#	this way avoids the Adjoint:
	y = fun(x', down[end])' # doing adjoint first
	y = fun(y, down[1])

#=	failed attempt at helping @inferred:
	f1 = (x, d) -> downsample1(x, d, warn=warn)
	tmp = x'
	tmp = hcat([f1(tmp[:,n], down[end]) for n=1:size(tmp,2)]...)'
	y = hcat([f1(tmp[:,n], down[1]) for n=1:size(tmp,2)]...)
@shows y
=#
=#

	return y
end

downsample2(x::AbstractArray{<:Number,2}, down::Int ; args...) =
	downsample2(x, [down, down] ; args...)


"""
    y = downsample3(x, down ; warn=true, T)

downsample by averaging by integer factors
in
- `x [nx ny nz]`
- `down` can be a scalar (same factor for all dimensions) or a 3-vector

option
- `warn::Bool` warn if noninteger multiple; default true
- `T::DataType` specify output eltype; default `eltype(x[1] / down[1])`

out
- `y [nx/down ny/down nz/down]`
"""
function downsample3(x::AbstractArray{<:Number,3},
		down::AbstractVector{Int}
		; warn::Bool = isinteractive(),
		T::DataType = eltype(x[1] / down[1])
	)

	length(down) != 3 && throw(DimensionMismatch("down $down"))

	idim = size(x)
	odim = floor.(Int, idim ./ down)

	warn && any(odim .* down .!= idim) && @warn("truncating to $odim")

	return downsample3_perm(x, Tuple(down)) # because it is faster
end

downsample3(x::AbstractArray{<:Number,3}, down::Int ; args...) =
	downsample3(x, [down, down, down] ; args...)


# this method is good for @inferred but is slower!
function downsample3_loop(x::AbstractArray{<:Number,3},
		down::AbstractVector{Int} ;
		T::DataType = eltype(x[1] / down[1]),
	)

	odim = floor.(Int, size(x) ./ down)

	y = Array{T}(undef, odim[1], odim[2], odim[3])
	d1 = down[1]
	d2 = down[2]
	d3 = down[3]
	d123 = d1 * d2 * d3

	for i3=1:odim[3]
		for i2=1:odim[2]
			for i1=1:odim[1]
				y[i1,i2,i3] = sum(@view x[(i1-1)*d1 .+ (1:d1),
					(i2-1)*d2 .+ (1:d2), (i3-1)*d3 .+ (1:d3)]) / d123
			end
		end
	end
	return y
end


# this method fails @inferred but is faster!
function downsample3_perm(x::AbstractArray{<:Number,3}, down::Dims{3})

	# down sample along each dimension
	y = downsample_dim1(x, down[1])
	y = downsample_dim1(permutedims(y, [2, 1, 3]), down[2])
	y = downsample_dim1(permutedims(y, [3, 2, 1]), down[3]) # [3 1 2] order
	y = permutedims(y, [2, 3, 1])

	return y
end
