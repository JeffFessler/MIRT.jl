#=
poweriter.jl
2019-06-06, Jeff Fessler, University of Michigan
=#

export poweriter

using LinearAlgebra: norm


"""
    v1,σ1 = poweriter(A; niter=?, ...)

Determine first right singular vector `v1`
and first singular value `σ1` of `A`
by applying power iteration to `A'A`

# in
- `A` M × N matrix

# option
- `niter` default 200
- `x0` initial guess of `v1`
- `tol` stopping tolerance for s1, default 1e-6
- `chat::Bool` verbose? default false

# out
- `v1` `[N]` principal right singular vector
- `σ1` spectral norm of `A`
"""
function poweriter(
    A ;
    niter::Int=200,
    tol::Real = 1e-6,
    x0::AbstractArray{<:Number} = ones(eltype(A), size(A,2)),
    chat::Bool = true,
)

    x = copy(x0)
    ratio_old = Inf
    for iter in 1:niter
        Ax = A * x
        ratio = norm(Ax) / norm(x)
        if abs(ratio - ratio_old) / ratio < tol
            chat && @info "done at iter $iter"
            break
        end
        ratio_old = ratio
        x = A' * Ax
        x /= norm(x)
    end
    return x, norm(A * x) / norm(x)
end
