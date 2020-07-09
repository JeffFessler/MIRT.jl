# poweriter.jl

using MIRT: poweriter

using Random: seed!
using LinearAlgebra: opnorm, norm
using Test: @test, @inferred


seed!(0)
M,N = 30,20

A = randn(M,N) # real
s0 = opnorm(A)
chat = false
_,s1 = @inferred poweriter(A ; tol=1e-9, chat=chat)
@test s0 ≈ s1

A = randn(ComplexF32, M, N) # complex
s0 = opnorm(A)
_,s1 = @inferred poweriter(A ; tol=1e-9, chat=chat)
@test s0 ≈ s1
