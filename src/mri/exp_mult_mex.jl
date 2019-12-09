using Random

include("max_percent_diff.jl")
export exp_mult_mex

"""
    D = exp_mult_mex(A, u, v)
in:
    A	[N L]	complex matrix
    u	[N]	vector
    v	[M]	vector
        one (and only one) of u and v must be complex!
out:
    D	[L M]	complex vector, D = A' * exp(-u * v.')
                D_lm = sum_n A_nl^* exp(-u_n v_m)
    
    This function is a memory efficient and fast implementation
    of AA matrix computation in int_tim_seg.m function in NUFFT package.
"""
function exp_mult_mex(A, u, v)
    A = A'
    return exp_mult_mex_helper(A, u, v)
end

function exp_mult_mex_helper(A, u::Vector{Complex{Float64}},  v::Vector{Float64})
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

# check this one, might not be correct for those
function exp_mult_mex_helper(A, u::Vector{Float64},  v::Vector{Complex{Float64}})
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

# check this one, might not be correct for those
function exp_mult_mex_helper(A, u::Vector{Float64},  v::Vector{Float64})
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

function exp_mult_mex_helper(A, u::Vector{Complex{Float64}},  v::Vector{Complex{Float64}})
    n = size(u)
    n = n[1]
    M = size(v)
    M = M[1]
    L, N = size(A)
    println("result returned, but not supposed have both u and v complex")

    D = zeros(Complex{Float64}, L, M)
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
    exp_mult_mex(test::Symbol)

Test function for exp_mult_mex.
Finds percent error of memory-efficient routine vs. normal routime
in
    test        Symbol  symbol to indicate test routine
out
    N/A
"""
function exp_mult_mex(test::Symbol)

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

    d1 = A' * exp.(-ur * transpose(v));
    d2 = exp_mult_mex(A, ur, v);
    println(max_percent_diff(d1, d2))

    d1 = A' * exp.(-u * transpose(vr));
    d2 = exp_mult_mex(A, u, vr);
    println(max_percent_diff(d1, d2))
end
