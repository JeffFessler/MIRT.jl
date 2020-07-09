# exp_mult.jl

using MIRT: exp_mult
using MIRT: exp_xform
using Test: @test


L = 10
N = 30
M = 20
A = randn(ComplexF32, N, L)
ur = randn(N)
ui = randn(N)
vr = randn(M)
vi = randn(M)
u = ur + im * ui
v = vr + im * vi

d1 = A' * exp.(-ur * transpose(vr))
@test d1 ≈ exp_mult(A, ur, vr)

d1 = A' * exp.(-ur * transpose(v))
@test d1 ≈ exp_mult(A, ur, v)

d1 = A' * exp.(-u * transpose(vr))
@test d1 ≈ exp_mult(A, u, vr)

d1 = A' * exp.(-u * transpose(v))
@test d1 ≈ exp_mult(A, u, v ; warnboth=false)

d1 = A' * exp.(-ur * transpose(v))
ur = reshape(ur, 1, :)
v = reshape(v, 1, :)
d2 = exp_xform(transpose(A'), ur, v, mode = :matrix)
d2 = transpose(d2)
@test d1 ≈ d2
