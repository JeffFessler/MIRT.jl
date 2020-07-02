#=
Afft.jl
2019-07-07, Jeff Fessler, University of Michigan
2020-06-30 in-place fft!
=#

export Afft

# using MIRT: embed!, getindex!
using LinearMapsAA: LinearMapAA, LinearMapAM, LinearMapAO
using FFTW: fft!, bfft!
using Test: @test


"""
    A = Afft(samp::AbstractArray{Bool} ; T::DataType = ComplexF32)
Make a LinearMapAO object for (under-sampled) FFT, of type `T`,
using given sampling pattern `samp`.
Especially for compressed sensing MRI with Cartesian sampling.

Option:
- `operator::Bool=true` set to `false` to return a `LinearMapAM`
- `work::AbstractArray` work space for in-place fft operations
"""
function Afft(samp::AbstractArray{Bool,D} ;
    T::DataType = ComplexF32,
    operator::Bool = true, # !
    work::AbstractArray{S,D} = Array{T}(undef, size(samp)),
) where {D, S <: Number}

    forw! = (y,x) -> getindex!(y, fft!(copyto!(work,x)), samp)
    back! = (x,y) -> bfft!(embed!(x,y,samp))

    dim = size(samp)

    if operator
        return LinearMapAA(forw!, back!,
            (sum(samp), prod(dim)), (name="fft",) ; T=T,
            idim = dim, # odim is 1D
            operator = operator,
        )
    end

    return LinearMapAA(
        (y,x) -> forw!(y, reshape(x,dim)),
    #   x -> fft(reshape(x,dim))[samp],
        (x,y) -> vec(back!(reshape(x,dim),y)),
    #   y -> vec(bfft(embed(y,samp))),
        (sum(samp), prod(dim)), (name="fft",) ; T=T,
    )
end


"""
    A = Afft(:test)
self test
"""
function Afft(test::Symbol)
    test != :test && throw(ArgumentError("test $test"))
    samp = trues(3,2); samp[2] = false

    @testset "AM" begin
        A = Afft(samp ; operator=false)
        @test A isa LinearMapAM
        @test Matrix(A)' ≈ Matrix(A')
        @test A * vec(ones(3,2)) == [6, 0, 0, 0, 0]
        @test A' * [1, 0, 0, 0, 0] == vec(ones(3,2))
        @test A.name == "fft"
        @test eltype(A) == ComplexF32
    end

    @testset "AO" begin
        A = Afft(samp)
        @test A isa LinearMapAO
        @test A * ones(3,2) == [6, 0, 0, 0, 0]
        @test A' * [1, 0, 0, 0, 0] == ones(3,2)
        @test Matrix(A)' ≈ Matrix(A')
    end

    true
end
