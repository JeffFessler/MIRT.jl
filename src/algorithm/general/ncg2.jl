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
    Tgf <: AbstractVector{<:Function},
    Tcf <: AbstractVector{<:Any},
    Tx <: AbstractArray{<:Number},
    Tp <: Any,
    Tu <: AbstractVector{<:AbstractArray},
    Tg <: AbstractArray{<:Number},
}
    B::Tb
    gradf::Tgf
    curvf::Tcf
    x::Tx # usually a Vector
    dir::Tx
    P::Tp
    Bx::Tu
    Bd::Tu
    grad_old::Tg
    grad_new::Tg
    npgrad::Tg
    betahow::Symbol
    iter::Int
    niter::Int
    ninner::Int

    function NCG(
        gradf::Tgf,
        curvf::Tcf,
        B::Tb,
        x::AbstractArray{<:Number,D}, # usually a Vector
        ;
        niter::Int = 10,
        ninner::Int = 5,
        P::Tp = I, # trick: this is an overloaded I (by LinearMapsAA)
        betahow::Symbol = :dai_yuan,
    ) where {
        Tb <: AbstractVector{<:Any},
        Tgf <: AbstractVector{<:Function},
        Tcf <: AbstractVector{<:Any},
        Tp <: Any,
        D
    }

        x = 1f0x # ensure >= Float32
        Tx = typeof(x)
        Bx = [Bj * x for Bj in B]
        Tu = typeof(Bx)

        # Gradients of each f_j must have same units.
        # The Tg line will fail (as it should) if units are incompatible.
        tmp = [gfj(Bxj) for (gfj, Bxj) in zip(gradf, Bx)]
        Tge = promote_type(eltype.(tmp)...) # cannot type infer tmp
# todo: force user to supply Tge?
        Tg = Array{Tge, D}

        axes(B) == axes(gradf) || axes(curvf) ||
            error("incompatible axes")
        new{Tb, Tgf, Tcf, Tx, Tp, Tu, Tg}(
            B, gradf, curvf, x, deepcopy(x), # dir
            P, Bx, deepcopy(Bx), # Bd
            similar(x, Tge), # grad_old
            similar(x, Tge), # grad_new
            similar(x, Tge), # npgrad
            betahow, 0, niter, ninner,
        )
    end
end


# Outer constructors


"""
    NCG(B, gradf, curvf, x0 ; out, ...)

Constructor for iterable
nonlinear preconditioned conjugate gradient algorithm
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
"""
function NCG(
    B::AbstractVector{<:Any},
    gradf::AbstractVector{<:Function},
    curvf::AbstractVector{<:Function},
    x0::AbstractArray{<:Number}, # usually a Vector
    ;
    kwargs...,
#   niter::Int = 50,
#   ninner::Int = 5,
#   P = I, # trick: this is an overloaded I (by LinearMapsAA)
#   betahow::Symbol = :dai_yuan,
)

    return NCG(gradf, curvf, B, x0; kwargs...)
end



"""
    ncg(args...; fun, out, kwargs...)

Convenience method for
nonlinear preconditioned conjugate gradient (NCG) algorithm.
See `NCG` documentation.

# option
- `fun(state)` User-defined function to be evaluated at each iteration.
- `out (niter+1) [fun(state), …, fun(state)]`
   * (all 0 by default). This is a Vector of length `niter+1`

# output
- `x` final iterate
"""

function ncg2(
    args...
#   B::AbstractVector{<:Any},
#   gradf::AbstractVector{<:Function},
#   curvf::AbstractVector{<:Function},
#   x0::AbstractArray{<:Number}, # usually a Vector
    ;
#   niter::Int = 50,
#   ninner::Int = 5,
#   P = I, # trick: this is an overloaded I (by LinearMapsAA)
#   betahow::Symbol = :dai_yuan,
    fun::Function = state -> 0, # todo: missing
    out::Union{Nothing,Vector{Any}} = nothing,
    kwargs...
)

    state = NCG(args... ; kwargs...)
    niter = state.niter
 
    out = Array{Any}(undef, niter+1)
    !isnothing(out) && length(out) < niter+1 && throw("length(out) < $niter+1")
    if !isnothing(out)
        out[1] = fun(state)
    end

    for item in state
        if !isnothing(out)
            out[state.iter+1] = fun(state)
        end
    end

    return state.x
end


function _update!(state::NCG)
    ninner = state.ninner
    B = state.B
    gradf = state.gradf
    curvf = state.curvf
    Bd = state.Bd
    Bx = state.Bx
    P = state.P
    x = state.x
    grad_old = state.grad_old
    grad_new = state.grad_new
    npgrad = state.npgrad
    dir = state.dir

    J = length(B)

    # todo: constructor
    grad = (Bx) -> sum([Bj' * gjf(Bjx) for (Bj, gjf, Bjx) in zip(B, gradf, Bx)])
    grad_new .= grad(Bx) # gradient: todo in place using Bd space

    mul!(npgrad, -P, grad_new)

    if state.iter == 0
        dir = npgrad
    else
        if state.betahow === :dai_yuan
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
    for (Bdj, Bj) in zip(Bd, B)
        mul!(Bdj, Bj, dir) # v_j in course notes
    end
    alf = line_search_mm(gradf, curvf, Bx, Bd; ninner) # todo: work
#   ls = LineSearchMM

    @. x += alf * dir # update x
    for (Bxj, Bdj) in zip(Bx, Bd) # update Bj * x
        @. Bxj += alf * Bdj
    end

    state.iter += 1

    return state
end

# Iterator

Base.IteratorSize(::NCG) = Base.SizeUnknown()
Base.IteratorEltype(::NCG) = Base.EltypeUnknown()

Base.iterate(state::NCG, arg=nothing) =
    (state.iter ≥ state.niter) ? nothing :
    (_update!(state), nothing)

Base.show(io::IO, ::MIME"text/plain", src::NCG) =
    _show_struct(io, MIME("text/plain"), src)
