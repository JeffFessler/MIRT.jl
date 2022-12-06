#=
Afft.jl
2019-07-07, Jeff Fessler, University of Michigan
2020-06-30 in-place fft!
2022-12-05 plan
=#

export Afft

# using MIRT: embed!, getindex!
using LinearMapsAA: LinearMapAA, LinearMapAM, LinearMapAO
using LinearAlgebra: mul!
using FFTW: plan_fft!, plan_bfft!


"""
    A = Afft(xdim::Dims, fdim ; T, operator, work, ...)
Make a `LinearMapAO` operator object
for FFT of size `xdim` along dimensions `fdim`, of type `T`.

In:
- `xdim::Dims{D}` dimensions of signal
- `fdim = 1:D` apply fft/bfft only along these dimensions

Option:
- `T::DataType = ComplexF32`
- `operator::Bool = true` set to `false` to return a `LinearMapAM`
- `work::AbstractArray` work space for in-place fft operations
- remaining arguments are passed to `plan_fft`

# Output
Returns a `LinearMapsAA.LinearMapA[M|O]` object.
"""
function Afft(
    xdim::Dims{D},
    fdim = 1:D,
    ;
    T::DataType = ComplexF32,
    operator::Bool = true, # !
    unitary::Bool = false,
    work::AbstractArray{S,D} = Array{T,D}(undef, xdim),
    kwargs...
) where {D}

    promote_type(S, T) == S || error("type S=$S cannot hold T=$T")
    size(work) == xdim || error("mismatch: work xdim")

    pf = plan_fft!(work, fdim; kwargs...)
    pb = plan_bfft!(work, fdim; kwargs...)
    if !fft_forward
        (pf, pb) = (pb, pf)
    end

    sqrtN = sqrt(prod(xdim))
    function forw!(y, x::AbstractArray{<:Real})
        @. work = complex(x)
        mul!(y, pf, work)
        unitary && y ./= sqrtN
        return y
    end
    function forw!(y, x)
        copyto!(y, x)
        mul!(y, pf, y)
        unitary && y ./= sqrtN
        return y
    end

    function back!(x, y)
        copyto!(x, y)
        unitary && x .*= sqrtN # todo merge copyto!
        mul!(x, pb, x)
        return x
    end

    prop = (; name="fft", fdim)

    if operator
        return LinearMapAA(forw!, back!,
            (1, 1) .* prod(dim), prop ; T,
            idim = xdim, odim = xdim,
            operator,
        )
    end

    return LinearMapAA(
        (y,x) -> forw!(y, reshape(x,xdim)), # fft(reshape(x,dim))
        (x,y) -> vec(back!(reshape(x,xdim),y)), # vec(bfft(y))
        (sum(samp), prod(dim)), prop ; T,
    )
end


"""
    A = Afft(samp::AbstractArray{Bool} ; T, dims, operator, work, ...)
Make a `LinearMapAO` operator object
for (under-sampled) FFT, of type `T`,
using given sampling pattern `samp`.
Especially for compressed sensing MRI with Cartesian sampling.

Option:
- `T::DataType = ComplexF32`
- `dims = 1:D` apply fft/bfft only along these dimensions
- `fft_forward::Bool = true` Use `false` to have `bfft!` in forward model.
- `operator::Bool = true` set to `false` to return a `LinearMapAM`
- `work::AbstractArray` work space for in-place fft operations
- remaining arguments are passed to `plan_fft`

# Output
Returns a `LinearMapsAA.LinearMapA[M|O]` object.
"""
function Afft(
    samp::AbstractArray{<:Bool, D},
    ;
    dims = 1:D,
    T::DataType = ComplexF32,
    operator::Bool = true, # !
    work::AbstractArray{S,D} = similar(samp, T),
    fft_forward::Bool = true,
    kwargs...
) where {D, S <: Number}

    promote_type(S, T) == S || error("type S=$S cannot hold T=$T")
    axes(work) == axes(samp) || error("axes mismatch: samp work")

    pf = plan_fft!(work, dims; kwargs...)
    pb = plan_bfft!(work, dims; kwargs...)
    if !fft_forward
        (pf, pb) = (pb, pf)
    end

#   forw! = (y,x) -> getindex!(y, fft!(copyto!(work,x)), samp)
    function forw!(y, x::AbstractArray{<:Real})
        @. work = complex(x)
        mul!(work, pf, work)
        getindex!(y, work, samp)
        return y
    end
    function forw!(y, x)
        copyto!(work, x)
        mul!(work, pf, work)
        getindex!(y, work, samp)
        return y
    end

#   back! = (x,y) -> bfft!(embed!(x,y,samp))
    function back!(x, y)
        embed!(x, y, samp)
        mul!(x, pb, x)
        return x
    end

    dim = size(samp)

    prop = (; name="fft", dims, fft_forward)

    if operator
        return LinearMapAA(forw!, back!,
            (sum(samp), prod(dim)), prop ; T,
            idim = dim, # odim is 1D
            operator,
        )
    end

    return LinearMapAA(
        (y,x) -> forw!(y, reshape(x,dim)), # fft(reshape(x,dim))[samp]
        (x,y) -> vec(back!(reshape(x,dim),y)), # vec(bfft(embed(y,samp)))
        (sum(samp), prod(dim)), prop ; T,
    )
end
