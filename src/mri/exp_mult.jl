#=
exp_mult.jl
translated from matlab/C file mri_exp_mult_mex.c
originally by Sangwoo Lee, University of Michigan, Jun, 2004
based on the 2003 file int_tim_seg.m by Brad Sutton
2004-6-21 modified by JF
2020-02 by Daniel Wan
=#

export exp_mult

using Test:@test

"""
    D = exp_mult(A, u, v)

Memory efficient and fast implementation of `D = A' * exp(-u * v^T)`
that is useful for B0-field-corrected MRI image reconstruction.

in:
* `A [N L]`	matrix
* `u [N]`	vector
* `v [M]`	vector
        usually exactly one of `u` and `v` are complex!

out:
* `D [L M]`	complex vector: `D = A' * exp(-u * v^T)`
`D_lm = sum_n A_nl^* exp(-u_n v_m)` 
"""
function exp_mult(A, u::AbstractVector{<:Number}, v::AbstractVector{<:Number})
    A = A'
    n = size(u,1)
    M = size(v,1)
    L, N = size(A)

    ~isreal(u) && ~isreal(v) && @warn "only one of u and v should be complex"

    T = promote_type(eltype(u), eltype(v), ComplexF32)
    D = zeros(T, L, M)
    col = Vector{T}(undef, N)
    
    for i = 1:M
        # compute ith column of N × M matrix B = exp(-u * v^T)
        for n = 1:N
            col[n] = exp(-u[n] * v[i])
        end
        # column of B computed
        for j = 1:L
            row = A[j, :]
            D[j,i] = transpose(row) * col
        end
    end
    return D
end


"""
    exp_mult(test::Symbol)

Test function for exp_mult.
Finds percent error of memory-efficient routine vs. normal routime
in
    test        Symbol  symbol to indicate test routine
out
    N/A
"""
function exp_mult(test::Symbol)
    L = 10
    N = 30
    M = 20
    A = randn(N,L) + im * randn(N,L)
    ur = randn(N)
    ui = randn(N)
    vr = randn(M)
    vi = randn(M)
    u = ur + im * ui
    v = vr + im * vi

    d1 = A' * exp.(-ur * transpose(vr))
    d2 = exp_mult(A, ur, vr);
    #println(isapprox(d1, d2))
    @test d1 ≈ d2

    d1 = A' * exp.(-ur * transpose(v))
    d2 = exp_mult(A, ur, v);
    #println(isapprox(d1, d2))
    @test d1 ≈ d2

    d1 = A' * exp.(-u * transpose(vr))
    d2 = exp_mult(A, u, vr);
    #println(isapprox(d1, d2))
    @test d1 ≈ d2

    d1 = A' * exp.(-u * transpose(v))
    d2 = exp_mult(A, u, v);
    #println(isapprox(d1, d2))
    @test d1 ≈ d2

    d1 = A' * exp.(-ur * transpose(v))
    ur = reshape(ur, 1, :)
    v = reshape(v, 1, :)
    d2 = exp_xform(transpose(A'), ur, v, mode = :matrix)
    d2 = transpose(d2)
    #println(isapprox(d1, d2))
    @test d1 ≈ d2

    true
end
