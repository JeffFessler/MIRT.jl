#=
Afft.jl
2019-07-07, Jeff Fessler, University of Michigan
2020-06-30 in-place fft!
=#

export Afft

# using MIRT: embed!, getindex!
using LinearMapsAA: LinearMapAA, LinearMapAM, LinearMapAO
using FFTW: fft!, bfft!


"""
    A = Afft(samp::AbstractArray{Bool} ; T::DataType = ComplexF32)
Make a LinearMapAO object for (under-sampled) FFT, of type `T`,
using given sampling pattern `samp`.
Especially for compressed sensing MRI with Cartesian sampling.

Option:
- `operator::Bool = true` set to `false` to return a `LinearMapAM`
- `work::AbstractArray` work space for in-place fft operations
"""
function Afft(
    samp::AbstractArray{Bool,D} ;
    T::DataType = ComplexF32,
    operator::Bool = true, # !
    work::AbstractArray{S,D} = similar(samp, T),
) where {D, S <: Number}

    forw! = (y,x) -> getindex!(y, fft!(copyto!(work,x)), samp)
    back! = (x,y) -> bfft!(embed!(x,y,samp))

    dim = size(samp)

    if operator
        return LinearMapAA(forw!, back!,
            (sum(samp), prod(dim)), (name="fft",) ; T,
            idim = dim, # odim is 1D
            operator,
        )
    end

    return LinearMapAA(
        (y,x) -> forw!(y, reshape(x,dim)), # fft(reshape(x,dim))[samp]
        (x,y) -> vec(back!(reshape(x,dim),y)), # vec(bfft(embed(y,samp)))
        (sum(samp), prod(dim)), (name="fft",) ; T,
    )
end
