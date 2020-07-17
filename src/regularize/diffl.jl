#=
diffl.jl
Left finite differences "in-place" (pre-allocated outpus)

Could use StaticKernels.jl for this; see timing test in ../../time.

Inspired by:
https://docs.julialang.org/en/latest/manual/performance-tips/#Pre-allocating-outputs-1

2019-06-22 Jeff Fessler, University of Michigan
=#

export diffl, diffl!, diffl_adj, diffl_adj!, diffl_map

using LinearMapsAA: LinearMapAA, LinearMapAM, LinearMapAO


"""
    diffl!(g::AbstractArray, x::AbstractArray, dim::Int ; ...)

Apply left finite difference operator to input array `x`,
storing the result "in-place" in pre-allocated output array `g`.

Arrays `g` and `x` must have the same size, and cannot alias.
The "first" elements of `g` are zero for dimension `dim`.
The default is `dim=1`.

Option:
- `add::Bool = false` use `x[i] + x[i-1]` instead of `x[i] - x[i-1]`
- `edge::Symbol = :zero` set the first elements of dimension `dim` to 0

Choose `edge=:circ` to use circulant (aka periodic) boundary conditions.
Choose `edge=:none` to leave the first elements untouched.

In 1D, if `x = [2, 6, 7]` then `g = [0, 4, 1]` with default options.
"""
function diffl!(
    g::AbstractArray{Tg,N},
    x::AbstractArray{Tx,N},
    dim::Int ;
    edge::Symbol=:zero,
    add::Bool=false,
) where {Tg,Tx,N}

    Base.require_one_based_indexing(g) && Base.require_one_based_indexing(x)
    1 <= dim <= N || throw(ArgumentError("dimension $dim out of range (1:$N)"))
    size(g) != size(x) && throw(DimensionMismatch("sizes g=>$(size(g)) vs x=>$(size(x))"))

    Nd = size(x, dim)
    if add
        @inbounds selectdim(g, dim, 2:Nd) .=
        selectdim(x, dim, 2:Nd) .+ selectdim(x, dim, 1:(Nd-1))
    else
        @inbounds selectdim(g, dim, 2:Nd) .=
        selectdim(x, dim, 2:Nd) .- selectdim(x, dim, 1:(Nd-1))
    end

    # handle edge conditions
    g1 = selectdim(g, dim, 1)
    if edge === :zero
        @inbounds g1 .= zero(Tg)
    elseif edge === :circ
        if add
            @inbounds g1 .= selectdim(x, dim, 1) .+ selectdim(x, dim, Nd)
        else
            @inbounds g1 .= selectdim(x, dim, 1) .- selectdim(x, dim, Nd)
        end
    else
        edge != :none && throw(ArgumentError("edge $edge"))
        # caution: in this case g1 is untouched, possibly undef
    end

    return g
end


# for default dim=1 case
diffl!(g::AbstractArray, x::AbstractArray ; kwargs...) = diffl!(g, x, 1 ; kwargs...)


"""
    diffl!(g::AbstractArray, x::AbstractArray, dims::AbstractVector{Int} ; ...)

When `x` is a `N`-dimensional array, the `i`th slice of the `g` array
(along its last dimension) is the `diffl!` of `x` along `dims[i]`.
This is useful for total variation (TV) and other regularizers
that need finite differences along multiple dimensions.
"""
function diffl!(
    g::AbstractArray{Tg,Ng},
    x::AbstractArray{Tx,Nx},
    dims::AbstractVector{Int} ;
    kwargs...,
) where {Tg,Tx,Ng,Nx}

    Ng != Nx+1 && throw(DimensionMismatch("Ng=$Ng Nx=$Nx"))
    size(g) != (size(x)..., length(dims)) &&
        throw(DimensionMismatch("sizes g=>$(size(g)) vs x=>$(size(x))"))

    for (i,d) in enumerate(dims)
        diffl!(selectdim(g, Ng, i), x, d ; kwargs...)
    end
    return g
end


"""
    g = diffl(x::AbstractArray ; ...)
Allocating version of `diffl!` along `dim=1`
"""
diffl(x::AbstractArray ; kwargs...) = diffl(x, 1 ; kwargs...)

"""
    g = diffl(x::AbstractArray, dim::Int ; ...)
Allocating version of `diffl!` along `dim`
"""
diffl(x::AbstractArray, dim::Int ; kwargs...) = diffl!(similar(x), x, dim ; kwargs...)

"""
    g = diffl(x::AbstractArray, dims::AbstractVector{Int} ; ...)
Allocating version of `diffl!` for `dims`
"""
diffl(x::AbstractArray, dims::AbstractVector{Int} ; kwargs...) =
    diffl!(similar(x, size(x)..., length(dims)), x, dims ; kwargs...)


"""
    diffl_adj!(z, g, dim::Int ; ...)

Adjoint of left finite difference `diffl!`, in-place.
Arrays `z` and `g` must be same size.
See `diffl` for details.
"""
function diffl_adj!(
    z::AbstractArray{Tz,N},
    g::AbstractArray{Tg,N},
    dim::Int ;
    reset0::Bool=true,
    edge::Symbol=:zero,
    add::Bool=false,
) where {Tz,Tg,N}

    1 <= dim <= N || throw(ArgumentError("dimension $dim out of range (1:$N)"))
    size(z) != size(g) && throw(DimensionMismatch("sizes z=>$(size(z)) vs g=>$(size(g))"))

    # todo: handle reset0 better
    if reset0
        z .= zero(Tz)
    end

    Nd = size(g, dim)

    if edge === :zero
        @inbounds selectdim(z, dim, 2:Nd) .+= selectdim(g, dim, 2:Nd)
    elseif edge === :circ
    #    selectdim(z, dim, 1:Nd) .+= selectdim(g, dim, 1:Nd)
        @inbounds z .+= g
    else
        edge != :none && throw(ArgumentError("edge $edge"))
        # in this case g1 is unused, even if undef
    end

    if add
        @inbounds selectdim(z, dim, 1:(Nd-1)) .+= selectdim(g, dim, 2:Nd)
        if edge === :circ
            @inbounds selectdim(z, dim, Nd) .+= selectdim(g, dim, 1)
        end
    else
        @inbounds selectdim(z, dim, 1:(Nd-1)) .-= selectdim(g, dim, 2:Nd)
        if edge === :circ
            @inbounds selectdim(z, dim, Nd) .-= selectdim(g, dim, 1)
        end
    end

    return z
end


"""
    diffl_adj!(z::AbstractArray, g::AbstractArray, dims::AbstractVector{Int} ; ...)

Adjoint of `diffl!` for multiple dimensions `dims`.
Here `g` must have one more dimension than `z`.
"""
function diffl_adj!(
    z::AbstractArray{Tz,Nz},
    g::AbstractArray{Tg,Ng},
    dims::AbstractVector{Int} ;
    kwargs...,
) where {Tz,Tg,Nz,Ng}

    Ng != Nz+1 && throw(DimensionMismatch("Ng=$Ng Nz=$Nz"))
    size(g) != (size(z)..., length(dims)) &&
        throw(DimensionMismatch("sizes g=>$(size(g)) vs z=>$(size(z))"))

    for (i,d) in enumerate(dims)
        diffl_adj!(z, selectdim(g, Ng, i), d ; reset0 = (i==1), kwargs...)
    end
    return z
end


"""
    z = diffl_adj(g::AbstractArray ; ...)
Allocating version of `diffl!` along `dim=1`
"""
diffl_adj(g::AbstractArray ; kwargs...) = diffl_adj(g, 1 ; kwargs...)

"""
    z = diffl(g::AbstractArray, dim::Int ; ...)
Allocating version of `diffl!` along `dim`
"""
diffl_adj(g::AbstractArray, dim::Int ; kwargs...) =
    diffl_adj!(similar(g), g, dim ; reset0=true, kwargs...)

"""
    z = diffl_adj(g::AbstractArray, dims::AbstractVector{Int} ; ...)
Allocating version of `diffl!` for `dims`
"""
function diffl_adj(
    g::AbstractArray{T,N},
    dims::AbstractVector{Int} ;
    kwargs...,
) where {T,N}

    size(g)[end] != length(dims) &&
        throw(DimensionMismatch("sizes g=>$(size(g)) vs dims=$dims"))
    return diffl_adj!(similar(g, size(g)[1:(N-1)]...), g, dims ; kwargs...)
end


"""
    T = diffl_map(N::Dims{D}, dims::AbstractVector{Int} ; kwargs...)
    T = diffl_map(N::Dims{D}, dim::Int ; kwargs...)

in
- `N::Dims` image size

options: see `diffl!`
- `T::Type` for `LinearMapAA`, default `Float32`
- `operator::Bool = true` use `false` for `LinearMapAM`

out
- `T` `LinearMapAA` object for computing finite differences via `T*x`
using `diffl!` and `diffl_adj!`
"""
function diffl_map(
    N::Dims{D},
    dims::AbstractVector{Int} ;
    T::Type=Float32,
    edge::Symbol = :zero,
    operator::Bool = true, # !
    kwargs...,
) where {D}

    !all(1 .<= dims .<= D) && throw(ArgumentError("dims $dims"))
    (edge == :none) && throw("edge=$edge unsupported")

    forw!(g,x) = diffl!(g, x, dims ; edge=edge, kwargs...)
    back!(z,g) = diffl_adj!(z, g, dims ; edge=edge, kwargs...)

    if operator
        return LinearMapAA(forw!, back!,
            (length(dims), 1) .* prod(N) ;
            prop = (name="diffl_map", N=N, dims=dims),
            T=T, operator=true,
            idim=N, odim=(N...,length(dims)),
        )
    end

    gshape = g -> reshape(g, N..., length(dims))
    return LinearMapAA(
        (g,x) -> forw!(gshape(g), reshape(x, N)),
        (z,g) -> vec(back!(reshape(z,N), gshape(g))),
        (length(dims), 1) .* prod(N),
        (name="diffl_map", N=N, dims=dims),
        T=T,
    )
end


# for single dimension case
function diffl_map(
    N::Dims{D},
    dim::Int ;
    T::Type = Float32,
    edge::Symbol = :zero,
    operator::Bool = true, # !
    kwargs...,
) where {D}

    (1 .<= dim .<= D) || throw(ArgumentError("dim $dim"))
    (edge == :none) && throw("edge=$edge unsupported")

    forw!(g,x) = diffl!(g, x, dim ; edge=edge, kwargs...)
    back!(z,g) = diffl_adj!(z, g, dim ; edge=edge, kwargs...)

    if operator
        return LinearMapAA(forw!, back!,
            (1, 1) .* prod(N) ;
            prop = (name="diffl_map", N=N, dim=dim),
            T=T, operator=true,
            idim=N, odim=N,
        )
    end

    return LinearMapAA(
        (g,x) -> forw!(reshape(g, N), reshape(x, N)),
        (z,g) -> vec(back!(reshape(z,N), reshape(g, N))),
        (1, 1) .* prod(N) ;
        prop = (name="diffl_map", N=N, dim=dim),
        T=T,
    )
end

# default dim=1
diffl_map(N::Dims ; kwargs...) = diffl_map(N, 1, ; kwargs...)
