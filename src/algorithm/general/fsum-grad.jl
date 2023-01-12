# fsum-grad.jl

using LinearAlgebra: mul!


"""
    _fsum_grad1!(grad, work, gradf, B, Bx)
Compute gradient of
``f(x) = \\sum_{j=1}^J f_j(B_j x)``
which is
``∇f(x) = \\sum_{j=1}^J B_j' ∇f_j(B_j x)``.

# in
- `grad` output gradient (mutated and returned)
- `work` workspace of same size as `x` (mutated)
- `gradf` ``∇f_1,…,∇f_J`` (single-argument gradient functions)
- `B` ``B_1,…,B_J``
- `Bx` ``B_1 x,…, B_J x``

For data with units,
both `grad` and `work` must have
`eltype` df/dx.

This function will be essentially non-allocating
if each ``∇f_j`` function has its own built-in workspace
or is otherwise memory efficient.
"""
function _fsum_grad1!(
    grad::AbstractArray{<:Number},
    work::AbstractArray{<:Number},
    gradf::AbstractVector{<:Function},
    B::AbstractVector{<:Any},
    Bx::AbstractVector{<:Any},
)
    J = length(gradf)
    mul!(grad, B[begin]', gradf[begin](Bx[begin]))
    for j in eachindex(gradf)[2:end]
        mul!(work, B[j]', gradf[j](Bx[j]))
        grad .+= work
    end
    return grad
end


# version for gradf!(grad, z)
# todo
function _fsum_grad2!(
    grad::AbstractArray{<:Number},
    work::AbstractArray{<:Number},
    gradf_work::AbstractVector{<:AbstractArray}, # ∇f_j(B_j x)
    gradf2::AbstractVector{<:Function},
    B::AbstractVector{<:Any},
    Bx::AbstractVector{<:Any},
)
    J = length(gradf)
    gradf2[begin](gradf_work[begin], Bx[begin])
    mul!(grad, B[begin]', gradf_work[begin])
    for j in eachindex(gradf)[2:end]
        gradf2[j](gradf_work[j], Bx[j])
        mul!(work, B[j]', gradf_work[j])
        grad .+= work
    end
    return grad
end


"""
    _fsum_make_grad(gradf, B)
Returns a function of `(grad, work, Bx)`
where `grad` and `work` are mutated.
See `_fsum_grad1!`.
"""
function _fsum_make_grad(
#   grad::AbstractArray{<:Number},
#   work::AbstractArray{<:Number},
    gradf::AbstractVector{<:Function},
    B::AbstractVector{<:Any},
)

    all(g -> _narg(g) == 1, gradf) || throw("not done")
    let B=B, gradf=gradf #, work=work, grad=grad
    #   (Bx) -> sum([B[j]' * gradf[j](Bx[j]) for j in 1:J])
        (grad, work, Bx) -> _fsum_grad1!(grad, work, gradf, B, Bx)
    end
end
