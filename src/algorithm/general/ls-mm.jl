#=
ls-mm.jl
Line-search based on majorize-minimize (MM) approach
=#

export line_search_mm

using LinearAlgebra: dot

_narg(fun::Function) = first(methods(fun)).nargs - 1

# This will fail (as it should) if the units of uj and vj are incompatible:
_ls_mm_worktype(uj, vj, α::Real) =
    typeof(oneunit(eltype(uj)) + α * oneunit(eltype(vj)))


"""
    LineSearchMMWork{Tz <: AbstractVector{<:AbstractArray}}

Workspace for storing ``z_j = u_j + α v_j`` in MM-based line search.

If all of those ``z_j`` arrays had the same `eltype`,
then we could save memory
by allocating just the longest vector needed.
But for Unitful data they could have different `eltype`s and `size`s,
which would require a lot of `reinterpret` and `reshape` to handle.
So we just allocate separate work arrays for each ``j``.
"""
mutable struct LineSearchMMWork{Tz <: AbstractVector{<:AbstractArray}}
    zz::Tz

    function LineSearchMMWork(
        uu::AbstractVector{<:AbstractArray},
        vv::AbstractVector{<:AbstractArray},
        α::Real,
    )

        axes(uu) == axes(vv) || error("incompatible u,v axes")
        all(j -> axes(uu[j]) == axes(vv[j]), eachindex(uu)) ||
            error("incompatible uj,vj axes")

        zz = [similar(uu[j], _ls_mm_worktype(uu[j], vv[j], α))
            for j in eachindex(uu)]
        Tz = typeof(zz)
        return new{Tz}(zz)
    end
end


"""
    LineSearchMM{...}

Mutable struct for MM-based line searches.
"""
mutable struct LineSearchMM{
    Tu <: AbstractVector{<:AbstractArray},
    Tv <: AbstractVector{<:AbstractArray},
    Tg <: AbstractVector{<:Function},
    Tc <: AbstractVector{<:Any},
    Tα <: Real,
    Tw <: LineSearchMMWork,
}
    uu::Tu # vector of B_j x
    vv::Tv # vector of B_j d
    dot_gradf::Tg # vector of <z, ∇f> functions
    dot_curvf::Tc # vector of <|z|², ω_f> functions
    α::Tα
    ninner::Int # max # of iterations
    iter::Int # initialized to 0
    work::Tw

    function LineSearchMM(
        uu::Tu,
        vv::Tv,
        dot_gradf::Tg,
        dot_curvf::Tc,
        α::Tα = 0f0,
        ninner::Int = 5,
        iter::Int = 0,
        work::Tw = LineSearchMMWork(uu, vv, α),
    ) where {
        Tu <: AbstractVector{<:AbstractArray},
        Tv <: AbstractVector{<:AbstractArray},
        Tg <: AbstractVector{<:Function},
        Tc <: AbstractVector{<:Any},
        Tα <: Real,
        Tw <: LineSearchMMWork,
    }

        all(==(axes(work.zz)), axes.((uu, vv, dot_gradf, dot_curvf))) ||
            error("incompatible axes")
        all(j -> axes(uu[j]) == axes(vv[j]), eachindex(uu)) ||
            error("incompatible u,v axes")
        all(j -> _ls_mm_worktype(uu[j], vv[j], α) == eltype(work.zz[j]),
            eachindex(uu)) || error("incompatible work type")

        Tαp = promote_type(Tα, typeof(1f0 * oneunit(Tα)))
        new{Tu, Tv, Tg, Tc, Tαp, Tw}(
            uu, vv, dot_gradf, dot_curvf, α, ninner, iter, work,
        )
    end
end


# Outer constructors


"""
    LineSearchMM(gradf, curvf, u, v; α0 ...)
    LineSearchMM(u, v, dot_gradf, dot_curvf; α0 ...)

Construct iterator for
line-search based on majorize-minimize (MM) approach
for a general family of 1D cost functions of the form
``h(α) = \\sum_{j=1}^J f_j(u_j + α v_j)``
where each function ``f_j(t)`` has a quadratic majorizer of the form
```math
q_j(t;s) = f_j(s) + ∇f_j(s) (t - s) + 1/2 ‖t - s‖^2_{C_j(s)}
```
where ``C_j(⋅)`` is diagonal matrix of curvatures.
(It suffices for each ``f_j`` to have a Lipschitz smooth gradient.)

Each function ``f_j : \\mathcal{X}_j ↦ \\mathbb{R}``
where conceptually
``\\mathcal{X}_j ⊆ \\mathbb{R}^{M_j}``,
but we allow more general domains.

There are two outer constructors (based on the positional arguments):
- The simple way (not type stable) provides
  - `gradf` vector of ``J`` functions return gradients of ``f_1,…,f_J``
  - `curvf` vector of ``J`` functions `z -> curv(z)` that return a scalar
    or a vector of curvature values for each element of ``z``
- The fancier way (type stable) provides
  - `dot_gradf::AbstractVector{<:Function} = make_dot_gradf.(gradf)`
    See `make_dot_gradf`.
  - `dot_curvf::AbstractVector{<:Function} = make_dot_curvf.(curvf)`
    See `make_dot_curvf`.

# in
- `u` vector of ``J`` arrays ``u_1,…,u_J`` (typically vectors)
- `v` vector of ``J`` arrays ``v_1,…,v_J`` (typically vectors)
We require `axes(u_j) == axes(v_j)` for all ``j=1,…,J``.

# option
- `α0::Real = 0f0` initial guess for step size
- `ninner::Int = 5` # max number of inner iterations of MM line search
- `work = LineSearchMMWork(u, v, α)` pre-allocated work space for ``u_j+α v_j``
"""
function LineSearchMM(
    uu::AbstractVector{<:AbstractArray},
    vv::AbstractVector{<:AbstractArray},
    dot_gradf::AbstractVector{<:Function},
    dot_curvf::AbstractVector{<:Function},
    ;
    α0::Real = 0f0,
    work::Tw = LineSearchMMWork(uu, vv, α0),
    ninner::Int = 5,
) where {Tw <: LineSearchMMWork}

    return LineSearchMM(uu, vv, dot_gradf, dot_curvf, α0, ninner, 0, work)
end


# This constructor is not type stable because of the broadcast:
function LineSearchMM(
    gradf::AbstractVector{<:Function},
    curvf::AbstractVector{<:Any},
    uu::AbstractVector{<:AbstractArray},
    vv::AbstractVector{<:AbstractArray},
    ; kwargs...
)
    dot_gradf = make_dot_gradf.(gradf)
    dot_curvf = make_dot_curvf.(curvf)

    return LineSearchMM(uu, vv, dot_gradf, dot_curvf; kwargs...)
end


#=
MM-based line search update for step size α
using h(α) = sum_j f_j(uj + α vj)
\dot{h}(α) = sum_j v_j' * ∇f_j(u_j + α v_j)
=#
function _update!(state::LineSearchMM)
    uu = state.uu
    vv = state.vv
    zz = state.work.zz
    α = state.α
    dot_gradf = state.dot_gradf
    dot_curvf = state.dot_curvf

    # Using Threads.@threads here slowed down the 3-ls-mm demo.
    for j in eachindex(zz)
        @. zz[j] = uu[j] + α * vv[j]
    end

    derh = sum(j -> real(dot_gradf[j](vv[j], zz[j])), eachindex(zz))
    curv = sum(j -> dot_curvf[j](vv[j], zz[j]), eachindex(zz))
    curv < zero(curv) && error("bug: curv=$curv < 0")
    if curv > zero(curv)
        state.α -= derh / curv
    end
    state.iter += 1
    return state
end


# Iterator

Base.IteratorSize(::LineSearchMM) = Base.SizeUnknown()
Base.IteratorEltype(::LineSearchMM) = Base.EltypeUnknown()

Base.iterate(state::LineSearchMM, arg=nothing) =
    (state.iter ≥ state.ninner) || (state.iter > 0 && iszero(state.α)) ? nothing :
    (_update!(state), nothing)


# Convenience wrapper

"""
    α = line_search_mm(args...; opt, fun, kwargs...)

Line-search based on majorize-minimize (MM) approach.
This is a wrapper around the iterator `LineSearchMM`.
See its constructors for `args` and other `kwargs`.

# option
- `fun(state)` User-defined function to be evaluated with the `state`
  initially and then after each iteration.
- `out::Union{Nothing,Vector{Any}} = nothing`
  optional place to store result of `fun` for iterates `0,…,ninner`:
   (All `missing by default.) This is a `Vector{Any}` of length `ninner+1`.

# output
- `α` final iterate

This function mutates the optional arguments `out` and `work`.
"""
function line_search_mm(
    args...
    ;
    out::Union{Nothing,Vector{Any}} = nothing,
    fun::Function = state -> missing,
    ninner::Int = 5,
    kwargs...
)

    !isnothing(out) && length(out) < ninner+1 && throw("length(out) < $(ninner+1)")

    state = LineSearchMM(args... ; ninner, kwargs...)

    if !isnothing(out)
        out[1] = fun(state)
    end

    for item in state
        if !isnothing(out)
             out[state.iter+1] = fun(state)
        end
    end

    return state.α
end


Base.show(io::IO, ::MIME"text/plain", src::LineSearchMM) =
    _show_struct(io, MIME("text/plain"), src)
Base.show(io::IO, ::MIME"text/plain", src::LineSearchMMWork) =
    _show_struct(io, MIME("text/plain"), src)
