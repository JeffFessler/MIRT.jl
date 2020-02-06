export exp_mult
using Test:@test
"""
    D = exp_mult(A, u, v)
in:
* `A	[N L]`	complex matrix
* `u	[N]`	vector
* `v	[M]`	vector
        one (and only one) of u and v must be complex!
out:
* `D	[L M]`	complex vector, D = A' * exp(-u * v.')
 `D_lm = sum_n A_nl^* exp(-u_n v_m)`
    
    This function is a memory efficient and fast implementation
    of AA matrix computation in int_tim_seg.m function in NUFFT package.
"""
function exp_mult(A, u, v)
    A = A'
    n = size(u)
    n = n[1]
    M = size(v)
    M = M[1]
    L, N = size(A)

    D = zeros(ComplexF64, L, M)
    col = Vector{Complex{Float64}}(undef, N)
    e = MathConstants.e
    for i = 1:M
        # compute ith column of B
        # B has M columns, N rows
        for j = 1:N
            col[j] = e ^ (-u[j] * v[i])
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
