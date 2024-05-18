#=
downsample.jl
Copyright 2019-03-05, Jeff Fessler, University of Michigan

see also MIRT/time/downsample.jl
=#

export downsample_dim1
export downsample1
export downsample2
export downsample3


"""
    y = downsample_dim1(x, down ; warn::Bool)

Down-sample `x` by factor `down` along first dimension by averaging.

# in
- `x [n1 (Nd)]`
- `down::Int` downsampling factor

# option
- `warn::Bool` warn if non-integer multiple; default `isinteractive()`

# out
- `y [n1÷down (Nd)]`
"""
function downsample_dim1(
    x::AbstractArray{<:Number},
    down::Int ;
    warn::Bool = isinteractive(),
)

    dim = size(x)
    dim1 = dim[1]
    x = reshape(x, dim1, :) # [n1 *Nd]
    m1 = dim1 ÷ down
    if m1 * down < dim1
        warn && @warn("truncating input size $dim1 to $(m1 * down)")
        x = x[1:(m1*down),:]
    end
    y = reshape(x, down, :)
    y = sum(y, dims=1) / down # mean
    y = reshape(y, m1, dim[2:end]...)
    return y
end


"""
    y = downsample1(x, down ; warn=true)

Downsample 1D vector by factor `down`.

# in
- `x [n1]`
- `down::Int` downsampling factor

# option
- `warn::Bool` warn if noninteger multiple; default `isinteractive()`

# out
- `y [n1/down]`
"""
function downsample1(
    x::AbstractVector{<:Number},
    down::Int ;
    warn::Bool = isinteractive(),
)

    dim = size(x)
    dim1 = dim[1]
    m1 = dim1 ÷ down
    if m1 * down < dim1
        warn && @warn("truncating input size $dim1 to $(m1 * down)")
        y = reshape((@view x[1:(m1*down)]), down, :)
    else
        y = reshape(x, down, :)
    end
    y = sum(y, dims=1) / down # mean [m1 * prod(dim[2:end])]
    return vec(y)
end


"""
    y = downsample2(x, down ; warn=true, T)

Downsample by averaging by integer factors.

# in
- `x [nx ny]`
- `down` can be a scalar (same factor for both dimensions) or a `NTuple{2,Int}`

# option
- `warn::Bool` warn if noninteger multiple; default `isinteractive()`
- `T::Type` specify output eltype; default `eltype(x[1] / down[1])`

# out
- `y [nx/down ny/down]`
"""
function downsample2(
    x::AbstractMatrix{<:Number},
    down::NTuple{2,Int} ;
    warn::Bool = isinteractive(),
    T::Type{<:Number} = eltype(x[1] / down[1])
)

    idim = size(x)
    odim = idim .÷ down

    warn && any(odim .* down .!= idim) && @warn("truncating to $odim")

    y = similar(x, T, odim)
    d1 = down[1]
    d2 = down[2]
    for i2 in 1:odim[2], i1 in 1:odim[1]
        y[i1,i2] =
            sum(@view x[(i1-1)*d1 .+ (1:d1), (i2-1)*d2 .+ (1:d2)]) / d1 / d2
    end

    fun = (x, d) -> downsample_dim1(x, d, warn=warn)

#=
#    this old way returns an adjoint type:
    y = fun(x, down[1])
    y = fun(y', down[2])'

#    this way avoids the adjoint:
    y = fun(x', down[2])' # doing adjoint first
    y = fun(y, down[1])
=#

    return y
end

downsample2(x::AbstractArray{<:Number,2}, down::Int ; args...) =
    downsample2(x, (down, down) ; args...)


"""
    y = downsample3(x, down ; warn=true, T)

Downsample by averaging by integer factors.

# in
- `x (nx,ny,nz)`
- `down` can be a scalar (same factor for all dimensions) or a `NTuple{3,Int}`

# option
- `warn::Bool` warn if noninteger multiple; default true
- `T::Type` specify output eltype; default `eltype(x[1] / down[1])`

# out
- `y (nx/down,ny/down,nz/down)`
"""
function downsample3(
    x::AbstractArray{<:Number,3},
    down::NTuple{3,Int} ;
    warn::Bool = isinteractive(),
    T::Type{<:Number} = eltype(x[1] / down[1]),
)

    idim = size(x)
    odim = idim .÷ down

    warn && any(odim .* down .!= idim) && @warn("truncating to $odim")

    return downsample3_perm(x, Tuple(down)) # because it is faster
end

downsample3(x::AbstractArray{<:Number,3}, down::Int ; args...) =
    downsample3(x, (down, down, down) ; args...)


# slower
function downsample3_loop(
    x::AbstractArray{<:Number,3},
    down::NTuple{3,Int} ;
    T::Type{<:Number} = eltype(x[1] / down[1]),
)

    odim = size(x) .÷ down

    y = similar(x, T, odim[1], odim[2], odim[3])
    d1 = down[1]
    d2 = down[2]
    d3 = down[3]
    d123 = d1 * d2 * d3

    for i3 in 1:odim[3], i2 in 1:odim[2], i1 in 1:odim[1]
        y[i1,i2,i3] = sum(@view x[(i1-1)*d1 .+ (1:d1),
            (i2-1)*d2 .+ (1:d2), (i3-1)*d3 .+ (1:d3)]) / d123
    end
    return y
end


# faster
function downsample3_perm(x::AbstractArray{<:Number,3}, down::NTuple{3,Int})

    # down sample along each dimension
    y = downsample_dim1(x, down[1])
    y = downsample_dim1(permutedims(y, [2, 1, 3]), down[2])
    y = downsample_dim1(permutedims(y, [3, 2, 1]), down[3]) # [3 1 2] order
    y = permutedims(y, [2, 3, 1])

    return y
end
