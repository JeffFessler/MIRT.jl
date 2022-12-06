#=
Asense.jl
=#

export Asense
using FFTW: fftshift!, ifftshift!, fft!, bfft!


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
- `fft_forward::Bool = true` Use `false` to have `bfft!` in forward model.

# Output
Returns a `LinearMapsAA.LinearMapAO` object.
"""
function Asense(
    samp::AbstractArray{<:Bool},
    smaps::Vector{<:AbstractArray{<:Number}},
    ;
    T::DataType = ComplexF32,
    fft_forward::Bool = true,
)

    Base.require_one_based_indexing(samp, smaps...)
    dims = size(samp)
    ncoil = length(smaps)
    for ic in 1:ncoil
        size(smaps[ic]) == dims || throw("size mismatch")
    end

    N = prod(dims)
    work1 = Array{T}(undef, dims)
    work2 = Array{T}(undef, dims)
    ffun! = fft_forward ? fft! : bfft!
    function forw!(y, x)
        for ic in 1:ncoil
            @. work1 = x * smaps[ic]
            fftshift!(work1, ffun!(ifftshift!(work2, work1)))
            y[:,ic] .= work1[samp]
        end
        return y
    end
    bfun! = fft_forward ? bfft! : fft!
    function back!(x, y)
        for ic in 1:ncoil
            embed!(work1, (@view y[:,ic]), samp)
            fftshift!(work1, bfun!(ifftshift!(work2, work1)))
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
    A = LinearMapAA(forw!, back!, (ncoil*count(samp), N);
        prop = (name = "Asense", samp, smaps),
        odim = (count(samp),ncoil), idim=dims, T,
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
