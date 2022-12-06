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
or a Vector of `ncoil` arrays of size `size(samp)`.

# Input
- `samp::AbstractArray{<:Bool}` `D`-dimensional sampling pattern.
- `smaps::Vector{<:AbstractArray{<:Number}}` or `AbstractArray{<:Number}`

# Option
- `T::DataType = ComplexF32`
- `dims = 1:D` apply fft/bfft only along these dimensions
- `fft_forward::Bool = true` Use `false` to have `bfft!` in forward model.
- remaining arguments are passed to `plan_fft`

# Output
Returns a `LinearMapsAA.LinearMapAO` object.
"""
function Asense(
    samp::AbstractArray{<:Bool, D},
    smaps::Vector{<:AbstractArray{<:Number}},
    ;
    dims = 1:D,
    T::DataType = ComplexF32,
    work1::AbstractArray{S,D} = similar(samp, T),
    work2::AbstractArray{S,D} = similar(samp, T),
    fft_forward::Bool = true,
    kwargs...
) where {D, S <: Number}

    promote_type(S, T) == S || error("type S=$S cannot hold T=$T")
    axes(work1) == axes(work2) == axes(samp) || error("axes mismatch: samp work")
    all(==(axes(samp)), axes.(smaps)) || throw("axes mismatch: samp smaps")

    sdim = size(samp)
    ncoil = length(smaps)

    pf = plan_fft!(work1, dims; kwargs...)
    pb = plan_bfft!(work1, dims; kwargs...)
    if !fft_forward
        (pf, pb) = (pb, pf)
    end

    function forw!(y, x)
        for ic in 1:ncoil
            @. work1 = x * smaps[ic]
#           fftshift!(work1, fft!(ifftshift!(work2, work1)))
            ifftshift!(work2, work1)
            mul!(work2, pf, work2)
            fftshift!(work1, work2)
            y[:,ic] .= work1[samp]
        end
        return y
    end

    function back!(x, y)
        for ic in 1:ncoil
            embed!(work1, (@view y[:,ic]), samp)
#           fftshift!(work1, bfft!(ifftshift!(work2, work1)))
            ifftshift!(work2, work1)
            mul!(work2, pb, work2)
            fftshift!(work1, work2)
            copyto!(work2, smaps[ic])
            conj!(work2)
            if ic == 1
                @. x = work1 * work2
            else
                @. x += work1 * work2
            end
        end
        return x
    end

    A = LinearMapAA(forw!, back!, (ncoil*count(samp), prod(sdim));
        prop = (name = "Asense", samp, smaps),
        odim = (count(samp),ncoil), idim=sdim, T,
    )
    return A
end


function Asense(
    samp::AbstractArray{<:Bool, D},
    smaps::AbstractArray{<:Number},
    ;
    kwargs...
) where D
    ndims(smaps) == D+1 || error("dimension mismatch")
    smapv = [eachslice(smaps, dims = D+1)...]
    return Asense(samp, smapv; kwargs...)
end
