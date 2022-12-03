#=
Asense.jl
=#

export Asense
using FFTW: fftshift!, ifftshift!, fft!, bfft!


"""
    Asense(samp, smaps; T)

Construct a MRI encoding matrix model
for Cartesian sampling pattern `samp`
and vector of coil sensitivity maps `smaps`.

Returns a `LinearMapAO` object.
"""
function Asense(
    samp::AbstractArray{<:Bool},
    smaps::Vector{<:Array{<:Number}},
    ;
    T::DataType = ComplexF32,
)

    Base.require_one_based_indexing(samp)
    dims = size(samp)
    ncoil = length(smaps)
    for ic in 1:ncoil
        size(smaps[ic]) == dims || throw("size mismatch")
    end

    N = prod(dims)
    work1 = Array{T}(undef, dims)
    work2 = Array{T}(undef, dims)
    function forw!(y, x)
        for ic in 1:ncoil
            @. work1 = x * smaps[ic]
            fftshift!(work1, fft!(ifftshift!(work2, work1)))
            y[:,ic] .= work1[samp]
        end
    end
    function back!(x, y)
        for ic in 1:ncoil
            embed!(work1, (@view y[:,ic]), samp)
            fftshift!(work1, bfft!(ifftshift!(work2, work1)))
            copyto!(work2, smaps[ic])
            conj!(work2)
            if ic == 1
                @. x = work1 * work2
            else
                @. x += work1 * work2
            end
        end
    end
    A = LinearMapAA(forw!, back!, (ncoil*count(samp), N);
        prop = (name = "Asense", samp, smaps),
        odim = (count(samp),ncoil), idim=dims, T,
    )
    return A
end
