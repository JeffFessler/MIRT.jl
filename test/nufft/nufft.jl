# test/nufft/nufft.jl

using MIRT: Anufft, nufft_init, nufft_errors
using MIRT: dtft_init
import MIRT: nufft_eltype

using NFFT
using LinearAlgebra: norm
using LinearMapsAA: LinearMapAM, LinearMapAO
using Test: @test, @testset, @test_throws


_norminf(a, b) = norm(b - a, Inf) / norm(a, Inf)


"""
    nufft_test1( ; M=30, N=20, n_shift=1.7, T=Float64, tol=7e-6)
Simple 1D NUFFT tests.
"""
function nufft_test1( ;
    M::Int = 30, N::Int = 20, n_shift::Real = 1.7,
    T::DataType = Float64, tol::Real = 7e-6,
)

    w = (rand(M) .- 0.5) * 2π
    w = T.(w)
    x = randn(complex(T), N)
    sd = dtft_init(w, N ; n_shift)
    (nufft, adjoint, A) = nufft_init(w, N ; n_shift, operator=false)
    o0 = sd.dtft(x)
    o1 = nufft(x)
    o2 = A * x
    y = Complex{T}.(o0 / norm(o0))
    b0 = sd.adjoint(y)
    b1 = adjoint(y)
    b2 = A' * y
    @test _norminf(o0, o1) < tol
    @test o2 == o1
    @test _norminf(b0, b1) < tol
    @test b2 == b1
    @test Matrix(A)' ≈ Matrix(A') # 1D adjoint test

    @test A.name == "nufft1"
    @test A isa LinearMapAM

    xint = ones(Int,N)
    i0 = sd.dtft(xint)
    i2 = A * xint
    @test _norminf(i0, i2) < tol

    # operation version
    (nufft, adjoint, B) = nufft_init(w, N; n_shift, do_many=false, operator=true)
    @test B isa LinearMapAO
    o3 = nufft(x)
    @test o3 == o1
    b3 = adjoint(y)
    @test b3 == b1

    # Anufft constructor
    A = Anufft(w, N; n_shift)
    @test A isa LinearMapAO
    @test A * x == o3
    @test A' * y == b3

    A = Anufft(w', N; n_shift) # w shape
    @test A._odim == (1,M)
    @test A * x == transpose(o3)
    @test A' * transpose(y) == b3

    nufft(ones(Int,N)) # produce a "conversion" warning
#   @test_throws MethodError nufft(ones(Int,N))
    true
end


"""
    nufft_test2( ; M=?, N=?, n_shift=?, T=?, tol=?)
Simple 2D NUFFT test.
"""
function nufft_test2( ;
    M::Int = 35,
    N::Dims = (10,8),
    n_shift::AbstractVector{<:Real} = [4,3],
    T::DataType = Float64, tol::Real = 2e-6,
)

    w = []

#=
    # fft sampling
    M = prod(N)
    w1 = 2π * (0:N[1]-1) / N[1]
    w2 = 2π * (0:N[2]-1) / N[2]
    w1 = repeat(w1, 1, N[2])
    w2 = repeat(w2', N[1], 1)
    w = [vec(w1) vec(w2)]
    n_shift = [0,0]
=#

    w = (rand(M,2) .- 0.5) * 2π

    w = T.(w)
    sd = dtft_init(w, N ; n_shift)
    (nufft, adjoint, A) = nufft_init(w, N ; n_shift, pi_error=false)

    x = randn(Complex{T}, N)
    o0 = sd.dtft(x)
    o1 = nufft(x)
    @test _norminf(o0, o1) < tol

#=
    # fft test only
    o2 = fft(x)
    o0 = reshape(o0, N)
    o1 = reshape(o1, N)
    @show _norminf(o0, o2)
    @show _norminf(o0, o1)
    @show _norminf(o1, o2)
=#

    y = convert(Array{Complex{T}}, o0 / norm(o0))
    b0 = sd.adjoint(y)
    b1 = adjoint(y)
    @test _norminf(b0, b1) < tol
    o2 = A * x
    @test o2 == o1
    b2 = A' * y
    @test b2 == b1
    @test Matrix(A)' ≈ Matrix(A') # 2D adjoint test

    @test nufft(cat(dims=4, x, 2x)) == cat(dims=3, nufft(x), nufft(2x))
    @test adjoint(cat(dims=3, y, 2y)) == cat(dims=4, adjoint(y), adjoint(2y))

    @test A.name == "nufft2"
    @test A.N == N

    Ao = nufft_init(w, N ; n_shift, pi_error=false, operator=true).A
    Am = nufft_init(w, N ; n_shift, pi_error=false, operator=false).A
    @test Ao * x == Am * vec(x)
    y = Ao * x
    @test Ao'*y == reshape(Am'*y, N)

    (nufft, adjoint, A) = nufft_init(w, N ; n_shift, pi_error=false,
        do_many=false, operator=true)
    o3 = nufft(x)
    @test _norminf(o0, o3) < tol

    A = Anufft(w, N ; n_shift)
    @test A isa LinearMapAO
    @test A * x == o3
# todo: adjoint

    if M == 35
        odim = (7, 5) # test array output
        w = reshape(w, odim..., 2)
        A = Anufft(w, N ; n_shift)
        @test A._odim == odim
        @test A * x == reshape(o3, odim)
    end

    true
end

@testset "errors" begin
    w, errs = nufft_errors( ; N=17)
    @test w isa AbstractRange
    @test errs isa Vector
end

@testset "basics" begin
    @test nufft_eltype(Bool) === Float32
    @test nufft_eltype(Float16) === Float32
    @test nufft_eltype(Float64) === Float64
    @test_throws String nufft_eltype(BigFloat)
    @test_throws ErrorException nufft_init([0], 2) # small
    @test_throws ArgumentError nufft_init([2π], 8) # π
end

@testset "1D" begin
    @test nufft_test1()
    @test nufft_test1(; N=41) # odd
    @test nufft_test1(; T=Float32)
end

@testset "2D" begin
    @test nufft_test2()
    @test nufft_test2(; T=Float32)
end


#=
    # todo: 1d vs 2d
    M = 4
    N = (M,1)
    w = (0:(M-1))/M * 2π
    w = [w zeros(M)]
    sd = dtft_init(w, N)
    Ad = Matrix(sd.A)
    sn1 = nufft_init(w[:,1], N[1]; pi_error=false)
    An1 = Matrix(sn1.A)
    @show maximum(abs.(An1 - Ad))
    sn2 = nufft_init(w, N; pi_error=false)
    An2 = Matrix(sn2.A)
#   @show o0 = sd.dtft(x)
#   @show o1 = sn.nufft(x)
#   @show o1-o0
#   @test maximum(abs.(o1 - o0)) < tol
=#


#=
    # todo MWE for 1D vs 2D
    M = 6
    x1 = collect((-(M÷2)):((M÷2)-1))/M
    N1 = M
    p1 = plan_nfft(x1, N1)
#   p1 = plan_nfft(x1', (N1,)) # this works too
#   f1 = ones(N1)
    f1 = [1.,2,3,4,0,0]
    o1 = nfft(p1, complex(f1)) # slightly annoying to have to use complex() here

    N2 = (M,4) # works fine with "4" instead of "1" here and separable signal
    x2 = [x1 zeros(M)]
    p2 = plan_nfft(x2', N2)
    f2 = f1 * [0,0,1,0]'
    o2 = nfft(p2, complex(f2)) # ditto

    display(round.(o1, digits=7))
    display(round.(o2, digits=7))

#   tmp = nufft_init(x1*2π, N1)
=#
