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
    A = Afft(xdim::Dims, fdim ; T, operator, unitary, work, ...)
Make a `<:LinearMapAX` operator object
for the FFT of an array of size `xdim`,
along dimensions `fdim`,
of type `T`.

In:
- `xdim::Dims{D}` array size
- `fdim = 1:D` apply fft/bfft only along these dimensions

Option:
- `T::DataType = ComplexF32`
- `operator::Bool = true` return a `LinearMapAO`
  set to `false` to return a `LinearMapAM`
- `unitary::Bool = false` set to `true` for unitary DFT
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
    work::AbstractArray{Tw,D} = Array{T,D}(undef, xdim),
    kwargs...
) where {D, Tw <: Number}

    all(in(1:D), fdim) || error("fdim $fdim")
    promote_type(Tw, T) == Tw || error("type Tw=$Tw cannot hold T=$T")
    size(work) == xdim || error("work size mismatch")

    # "work" is essential below because the plans are designed
    # for eltype T only and x,y might have different eltype.
    pf = plan_fft!(work, fdim; kwargs...)
    pb = plan_bfft!(work, fdim; kwargs...)

    N = prod(xdim[fdim])

    function fun!(y::AbstractArray{Ty}, x, plan, factor::Real,
    ) where {Ty <: Number}
        if Ty == Tw
            if factor != 1
                @. y = x * factor
            else
                copyto!(y, x)
            end
            mul!(y, plan, y) # no "work" needed in this case
        else
            if factor != 1
                @. work = x * factor
            else
                copyto!(work, x)
            end
            mul!(work, plan, work)
            copyto!(y, work)
        end
        return y
    end

    Tr = real(T)
    factor = unitary ? Tr(1/sqrt(N)) : 1
    forw!(y, x) = fun!(y, x, pf, factor)
    back!(x, y) = fun!(x, y, pb, factor)

    prop = (; name="fft", fdim, unitary, pf, pb)

    if operator
        return LinearMapAA(forw!, back!,
            (1, 1) .* prod(xdim), prop ; T,
            idim = xdim, odim = xdim,
            operator,
        )
    end

    return LinearMapAA(
        (y,x) -> vec(forw!(reshape(y,xdim), reshape(x,xdim))),
        (x,y) -> vec(back!(reshape(x,xdim), reshape(y,xdim))),
        (1, 1) .* prod(xdim), prop ; T,
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
    work::AbstractArray{Tw,D} = similar(samp, T),
    fft_forward::Bool = true,
    kwargs...
) where {D, Tw <: Number}

    promote_type(Tw, T) == Tw || error("type Tw=$Tw cannot hold T=$T")
    axes(work) == axes(samp) || error("axes mismatch: samp work")

    pf = plan_fft!(work, dims; kwargs...)
    pb = plan_bfft!(work, dims; kwargs...)
    if !fft_forward
        (pf, pb) = (pb, pf)
    end

#   forw! = (y,x) -> getindex!(y, fft!(copyto!(work,x)), samp)
    function forw!(y, x)
        copyto!(work, x)
        mul!(work, pf, work)
        getindex!(y, work, samp)
        return y
    end

#   back! = (x,y) -> bfft!(embed!(x,y,samp))
    function back!(x::AbstractArray{<:Tx}, y) where {Tx <: Number}
        if Tx == Tw
            embed!(x, y, samp)
            mul!(x, pb, x) # no work needed
        else
            embed!(work, y, samp)
            mul!(work, pb, work)
            copyto!(x, work)
        end
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
