#=
errors.jl
Worst-case errors for NUFFT
=#

export nufft_errors

#using MIRT: dtft_init, nufft_init
#include("dtft.jl")
using NFFT: plan_nfft, nfft, nfft_adjoint
using LinearAlgebra: norm



"""
    w, errs = nufft_errors( ; M=401, w=?, N=513, n_shift=0, ...)

Compute NUFFT approximation errors (for signal of length `N` of unit norm),
for given digital frequency values `w`, i.e., Ω.
Default `w` is `range(0, 2π/N, M)`.
"""
function nufft_errors( ;
    M::Int = 401,
    N::Int = 512,
    w::AbstractArray{<:Real} = LinRange(0, 2π/N, M), # Ω values
    n_shift::Real = 0,
    kwargs...,
)

    sd = dtft_init(w, N ; n_shift)
    sn = nufft_init(w, N ; n_shift, kwargs...)
    E = Matrix(sn.A) - Matrix(sd.A) # (M,N)
    return w, vec(mapslices(norm, E, dims=2)) # [M]
end
