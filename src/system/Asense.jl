#=
Asense.jl
=#

export Asense

# using MIRT: embed!, getindex!
using FFTW: fftshift!, ifftshift!, fft!, bfft!
using LinearAlgebra: mul!
using FFTW: plan_fft!, plan_bfft!


"""
    Asense(samp, smaps; T)

Construct a MRI encoding operator model
for `D`-dimensional Cartesian sampling pattern `samp`
and coil sensitivity maps `smaps`.

The input `smaps` can either be a `D+1` dimensional array
of size `(size(samp)..., ncoil)`,
or a Vector (or `Slices`) of `ncoil` arrays of size `size(samp)`.

# Input
- `samp::AbstractArray{<:Bool}` `D`-dimensional sampling pattern.
- `smaps::Vector{<:AbstractArray{<:Number}}` or `AbstractArray{<:Number}`

# Option
- `T::Type = ComplexF32`
- `dims = 1:D` apply fft/bfft only along these dimensions
- `fft_forward::Bool = true` Use `false` to have `bfft!` in forward model.
- `unitary::Bool = false` set to `true` for unitary DFT
- remaining arguments are passed to `plan_fft`

# Output
Returns a `LinearMapsAA.LinearMapAO` object.
"""
function Asense(
    samp::AbstractArray{<:Bool, D},
    smaps::AbstractVector{<:AbstractArray{<:Number}},
    ;
    dims = 1:D,
    T::Type{<:Complex{<:AbstractFloat}} = ComplexF32,
    work1::AbstractArray{Tw,D} = similar(samp, T),
    work2::AbstractArray{Tw,D} = similar(samp, T),
    unitary::Bool = false,
    fft_forward::Bool = true,
    kwargs...
) where {D, Tw <: Number}

    all(in(1:D), dims) || throw(DimensionMismath("dims $dims"))
    promote_type(Tw, T) == Tw || throw("type Tw=$Tw cannot hold T=$T")
    axes(work1) == axes(work2) == axes(samp) || throw("axes mismatch: samp work")
    all(==(axes(samp)), axes.(smaps)) || throw("axes mismatch: samp smaps")

    sdim = size(samp)
    ncoil = length(smaps)

    pf = plan_fft!(work1, dims; kwargs...)
    pb = plan_bfft!(work1, dims; kwargs...)
    if !fft_forward
        (pf, pb) = (pb, pf)
    end

    factor = unitary ? 1/sqrt(prod(size(samp)[dims])) : 1

    function forw!(y, x)
        for ic in 1:ncoil
            @. work1 = x * smaps[ic] # apply sensitivity map
            ifftshift!(work2, work1)
            mul!(work2, pf, work2) # FFT
            fftshift!(work1, work2)
            if factor == 1
                y[:,ic] .= work1[samp] # sampling
            else
                @. y[:,ic] = work1[samp] * factor
            end
        end
        return y
    end

    function back!(x, y)
        for ic in 1:ncoil
            embed!(work1, (@view y[:,ic]), samp)
            ifftshift!(work2, work1)
            mul!(work2, pb, work2)
            fftshift!(work1, work2)
            copyto!(work2, smaps[ic])
            if ic == 1
                if factor == 1
                    @. x = work1 * conj(work2)
                else
                    @. x = work1 * conj(work2) * factor
                end
            else
                if factor == 1
                    @. x += work1 * conj(work2)
                else
                    @. x += work1 * conj(work2) * factor
                end
            end
        end
        return x
    end

    A = LinearMapAA(forw!, back!, (ncoil*count(samp), prod(sdim));
        prop = (name = "Asense", samp, smaps, unitary, fft_forward),
        odim = (count(samp),ncoil), idim=sdim, T,
    )
    return A
end


# handle typical array version of `smaps` efficiently
function Asense(
    samp::AbstractArray{<:Bool, D},
    smaps::AbstractArray{<:Number},
    ;
    kwargs...
) where D
    ndims(smaps) == D+1 || throw(DimensionMismatch("$(ndims(smaps)) ≠ $(D+1)"))
    smapv = eachslice(smaps, dims = D+1) # Slices
    return Asense(samp, smapv; kwargs...)
end
