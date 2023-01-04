#=
ncg2.jl
Nonlinear CG optimization (mutating version)
2023-01-03, Jeff Fessler, University of Michigan
=#

export ncg2

using LinearAlgebra: I, dot
# using MIRT: line_search_mm, _show_struct


"""
    NCG{...}

Mutable struct for nonlinear CG.
"""
mutable struct NCG{
    Tb <: AbstractVector{<:Any},
    Tg <: AbstractVector{<:Function},
    Tc <: AbstractVector{<:Any},
    Tx <: AbstractArray{<:Number},
    Tp <: Any,
    Tu <: AbstractVector{<:Any},
}
    B::Tb
    gradf::Tg
    curvf::Tc
    x::Tx # usually a Vector
    P::Tp
    Bx::Tu
    Bd::Tu
    niter::Int
#   ninner::Int

    function NCG(
        B::Tb,
        gradf::Tg
        curvf::Tc,
        x0::Tx # usually a Vector
        ;
        niter::Int = 50,
#       ninner::Int = 5,
        P::Tp = I, # trick: this is an overloaded I (by LinearMapsAA)
        betahow::Symbol = :dai_yuan,
    ) where {
        Tb <: AbstractVector{<:Any},
        Tg <: AbstractVector{<:Function},
        Tc <: AbstractVector{<:Any},
        Tx <: AbstractArray{<:Number},
        Tp <: Any,
    }

        Bx = [B[j] * x for j in 1:J]

        axes(B) == axes(gradf) || axes(curvf) ||
            error("incompatible axes")
        new{Tb, Tg, Tc, Tx, Tp, Tu}(
            B, gradf, curvf, x0, P, Bx, deepcopy(Bx), niter,
        )

    end

end


# Outer constructors


"""
    (x,out) = ncg(B, gradf, curvf, x0 ; ...)

Nonlinear preconditioned conjugate gradient algorithm
to minimize a general "inverse problem" cost function of the form
``\\Psi(x) = \\sum_{j=1}^J f_j(B_j x)``
where each function ``f_j(t)`` has a quadratic majorizer of the form
```math
q_j(t;s) = f_j(s) + ∇f_j(s) (t - s) + 1/2 ‖t - s‖^2_{C_j(s)}
```
where ``C_j(⋅)`` is diagonal matrix of curvatures.
(It suffices for each ``f_j`` to have a Lipschitz smooth gradient.)

This CG method uses a majorize-minimize (MM) line search.

# in
- `B` vector of ``J`` blocks ``B_1,…,B_J``
- `gradf` vector of ``J`` functions return gradients of ``f_1,…,f_J``
- `curvf` vector of ``J`` functions `z -> curv(z)` that return a scalar
  or a vector of curvature values for each element of ``z``
- `x0` initial guess; need `length(x) == size(B[j],2)` for ``j=1,…,J``

Usually `x0` is a `Vector` but it can be an `Array`
if each `B_j` is a linear operator (e.g., `LinearMapAO`)
of suitable "dimensions".

# option
- `niter = 50` # number of outer iterations
- `ninner = 5` # number of inner iterations of MM line search
- `P = I` # preconditioner
- `betahow = :dai_yuan` "beta" method for the search direction
- `fun` User-defined function to be evaluated with two arguments (x,iter).
   * It is evaluated at `(x0,0)` and then after each iteration.

# output
- `x` final iterate
- `out` `[niter+1] (fun(x0,0), fun(x1,1), ..., fun(x_niter,niter))`
   * (all 0 by default). This is an array of length `niter+1`
"""
function ncg2(
    B::AbstractVector{<:Any},
    gradf::AbstractVector{<:Function},
    curvf::AbstractVector{<:Function},
    x0::AbstractArray{<:Number}, # usually a Vector
    ;
#   out::Union{Nothing,Vector{Any}} = nothing,
    niter::Int = 50,
    ninner::Int = 5,
    P = I, # trick: this is an overloaded I (by LinearMapsAA)
    betahow::Symbol = :dai_yuan,
    fun::Function = (x,iter) -> 0, # todo: missing
)

    Base.require_one_based_indexing(B, gradf, curvf)

    out = Array{Any}(undef, niter+1)
    !isnothing(out) && length(out) < niter && throw("length(out) < $niter")
    if !isnothing(out)
        out[1] = fun(x0, 0) # todo: state
    end

    for item in state
    end

    return state.x #, out (todo)
#   return eltype(x0).(x), out # todo
end


function _update!(state::NCG)
    B = state.B
    Bd = state.Bd
    Bx = state.Bx
    grad_old = state.grad_old
    grad_new = state.grad_new
    npgrad = state.npgrad
    dir = state.dir

    J = length(B)

    dir = []
    grad_new = []

    grad = (Bx) -> sum([B[j]' * gradf[j](Bx[j]) for j in 1:J]) # todo: constructor

    grad_new .= grad(Bx) # gradient: todo in place
    mul!(npgrad, -P, grad_new)
    if state.iter == 0
        dir = npgrad
    else
        if betahow === :dai_yuan
            denom = dot(grad_new, dir) - dot(grad_old, dir)
            if iszero(denom)
                betaval = 0
            else
#               betaval = dot(grad_new, P * grad_new) / denom
                betaval = -dot(grad_new, npgrad) / denom
            end
        else
            throw(ArgumentError("unknown beta choice: $betahow"))
        end
        @. dir = npgrad + betaval * dir # search direction
    end
    grad_old .= grad_new

    # MM-based line search for step size alpha
    # using h(a) = sum_j f_j(uj + a vj)
    for j in 1:J
        mul!(Bd[j], @view B[j], dir) # v_j in course notes, todo: view?
    end
    alf = line_search_mm(gradf, curvf, Bx, Bd; state.ninner) # todo: work
#   ls = LineSearchMM(

    @. x += alf * dir # update x
    for j in 1:J # update Bj * x
        @. Bx[j] += alf * Bd[j]
    end

    return state
end

# Iterator

Base.IteratorSize(::NCG) = Base.SizeUnknown()
Base.IteratorEltype(::NCG) = Base.EltypeUnknown()

Base.iterate(state::NCG, arg=nothing) =
    (state.iter ≥ state.niiter) ? nothing :
    (_update!(state), nothing)

Base.show(io::IO, ::MIME"text/plain", src::NCG) =
    _show_struct(io, MIME("text/plain"), src)
