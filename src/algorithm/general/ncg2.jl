#=
ncg2.jl
Nonlinear CG optimization (mutating version)
2023-01-03, Jeff Fessler, University of Michigan
=#

export ncg2

using LinearAlgebra: I, dot
# using MIRT: line_search_mm, LineSearchMMWork, _show_struct

Base.allequal(f::Function, itr) = allequal(f.(itr))


"""
    NCG{...}

Mutable struct for nonlinear CG.
"""
mutable struct NCG{
    Tb <: AbstractVector{<:Any},
    Tgf <: AbstractVector{<:Function},
    Tdg <: AbstractVector{<:Function}, # dot_gradf
    Tdc <: AbstractVector{<:Function}, # dot_curvf
    Tg! <: Function, # grad!
    Tx <: AbstractArray{<:Number},
    Tp <: Any,
    Tu <: AbstractVector{<:AbstractArray},
    Tg <: AbstractArray{<:Number},
#   Tl <: LineSearchMMWork,
    Tl <: LineSearchMM,
}
    B::Tb
    gradf::Tgf
    dot_gradf::Tdg
    dot_curvf::Tdc
    grad!::Tg!
    x::Tx # usually a Vector
    dir::Tx
    P::Tp
    Bx::Tu
    Bd::Tu
    grad_old::Tg
    grad_new::Tg
    npgrad::Tg
#   ls_work::Tl
    lsmm::Tl
    betahow::Symbol
    iter::Int
    niter::Int

    """
todo
    """
    function NCG(
        B::Tb,
        gradf::Tgf,
        curvf::Tcf, # used only for dot_curvf
        x::AbstractArray{<:Number,D}, # usually a Vector
        ;
        niter::Int = 10,
        ninner::Int = 5,
        P::Tp = I, # trick: this is an overloaded I (by LinearMapsAA)
        betahow::Symbol = :dai_yuan,
        dot_gradf::Tdg = make_dot_gradf.(gradf),
        dot_curvf::Tdc = make_dot_curvf.(curvf),
    ) where {
        Tb <: AbstractVector{<:Any},
        Tgf <: AbstractVector{<:Function},
        Tcf <: AbstractVector{<:Any}, # think Union{Function,RealU}
        Tdg <: AbstractVector{<:Function},
        Tdc <: AbstractVector{<:Function},
        Tp <: Any,
        D
    }

        allequal(axes, [gradf, curvf, dot_gradf, dot_curvf]) ||
            error("incompatible axes")

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

        Bd = deepcopy(Bx)
        lsmm = LineSearchMM(Bx, Bd, dot_gradf, dot_curvf; ninner)
        # @show Base.summarysize(lsmm) / 1024^2
        # lsmm.work also allocates something the size of Bx
        # so all together there are 3 such arrays.

        dir = similar(x)
        npgrad = similar(x, Tge)
        # trick: use "dir" as space for grad, npgrad as work
#       grad! = _fsum_make_grad(dir, npgrad, gradf, B)
#grad! = _fsum_make_grad(similar(x), similar(x), gradf, B)
grad! = _fsum_make_grad(gradf, B)
        Tg! = typeof(grad!)

#       ls_work = LineSearchMMWork(Bx, Bx, 0f0)
#       Tl = typeof(ls_work)
        Tl = typeof(lsmm)

        new{Tb, Tgf, Tdg, Tdc, Tg!, Tx, Tp, Tu, Tg, Tl}(
            B, gradf, dot_gradf, dot_curvf,
            grad!, x, dir, P, Bx, Bd,
            similar(x, Tge), # grad_old
            similar(x, Tge), # grad_new # todo: use with grad?
            npgrad,
#           ls_work,
            lsmm,
            betahow, 0, niter,
        )
    end
end


# Outer constructors

tmp =
"
    NCG(B, gradf, curvf, x0 ; out, ...)

Constructor for iterable
nonlinear preconditioned conjugate gradient algorithm
to minimize a general 'inverse problem' cost function of the form
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
of suitable 'dimensions'.

# option
- `niter = 50` # number of outer iterations
- `ninner = 5` # number of inner iterations of MM line search
- `P = I` # preconditioner
- `betahow = :dai_yuan` 'beta' method for the search direction
"

@doc tmp NCG

#=
function NCG(
    B::AbstractVector{<:Any},
    gradf::AbstractVector{<:Function},
    curvf::AbstractVector{<:Any},
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
=#


# Mutating update


function _update!(state::NCG)
    B = state.B
    Bd = state.Bd
    Bx = state.Bx
    x = state.x
    grad_old = state.grad_old
    grad_new = state.grad_new
    npgrad = state.npgrad
    dir = state.dir

    state.grad!(grad_new, npgrad, Bx) # trick: use npgrad as work

    if state.P === I
        @. npgrad = grad_new
    else
        mul!(npgrad, state.P, grad_new)
        @. npgrad = -npgrad # -P * grad_new
    end

    if state.iter == 0
        dir .= npgrad
    else
        if state.betahow === :zero # GD
            betaval = 0
        elseif state.betahow === :dai_yuan
            denom = dot(grad_new, dir) - dot(grad_old, dir)
            if iszero(denom)
                betaval = 0
            else
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
        mul!(Bdj, Bj, dir) # (B_j * dir)
    end

    state.lsmm.iter = 0
    state.lsmm.α = 0 # uu=Bx and vv=Bd already set up!
    for _ in state.lsmm # run LS-MM
    end
    alf = state.lsmm.α
#=
    alf = line_search_mm(
        Bx, Bd,
        state.dot_gradf, state.dot_curvf;
        state.ninner, work = state.ls_work,
    )
=#

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

Base.iterate(state::NCG, arg = nothing) =
    (state.iter ≥ state.niter) ? nothing :
    (_update!(state), nothing)


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
function ncg2( # todo
    args...
    ;
    fun::Function = isnothing,
    out::Union{Nothing,Vector{Any}} = nothing,
    kwargs...
)

    state = NCG(args... ; kwargs...)
    niter = state.niter

    !isnothing(out) && length(out) < niter+1 && throw("length(out) < $niter+1")
    !isnothing(out) && fun == isnothing && throw("`out` requires `fun`")
    if !isnothing(out)
        out[1] = fun(state)
    end

    for _ in state
        if !isnothing(out)
            out[state.iter+1] = fun(state)
        end
    end

    return state.x
end

Base.show(io::IO, ::MIME"text/plain", src::NCG) =
    _show_struct(io, MIME("text/plain"), src)
