#=
exp_mult.jl
translated from matlab/C file mri_exp_mult_mex.c
originally by Sangwoo Lee, University of Michigan, Jun, 2004
based on the 2003 file int_tim_seg.m by Brad Sutton
2004-6-21 modified by JF
2020-02 by Daniel Wan
=#

export exp_mult


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
