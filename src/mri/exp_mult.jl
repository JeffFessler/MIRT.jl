#=
exp_mult.jl
translated from matlab/C file mri_exp_mult_mex.c
originally by Sangwoo Lee, University of Michigan, Jun, 2004
based on the 2003 file int_tim_seg.m by Brad Sutton
2004-6-21 modified by JF
2020-02 by Daniel Wan
=#

export exp_mult

#using MIRT: exp_xform
using Test:@test


"""
    D = exp_mult(A, u, v ; warnboth)

Memory efficient and fast implementation of `D = A' * exp(-u * v^T)`
that is useful for B0-field-corrected MRI image reconstruction.

in:
* `A [N L]` matrix
* `u [N]` vector
* `v [M]` vector
* `warnboth` warn if both `u` and `v` are complex; default: true

out:
* `D [L M]` complex vector: `D = A' * exp(-u * v^T)`
`D_lm = sum_n A_nl^* exp(-u_n v_m)`
"""
function exp_mult(A, u::AbstractVector{<:Number}, v::AbstractVector{<:Number}
        ; warnboth::Bool = true)

    N, L = size(A)
    n = size(u,1)
    M = size(v,1)

    ~isreal(u) && ~isreal(v) && warnboth &&
        @warn "only one of u and v should be complex"

    T = promote_type(eltype(u), eltype(v), ComplexF32)
    D = zeros(T, L, M)

    for m = 1:M
        col = exp.(-u * v[m]) # [N] mth column of N × M matrix B = exp(-u * v^T)
        D[:,m] = A' * col
    end
    return D
end


"""
    exp_mult(:test)
self test
"""
function exp_mult(test::Symbol)
    test != :test && throw("bad test $test")

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

    true
end
