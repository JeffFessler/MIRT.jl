#=
ls-mm.jl
Line-search based on majorize-minimize (MM) approach
=#

export line_search_mm

using LinearAlgebra: dot

_narg(fun::Function) = first(methods(fun)).nargs - 1


"""
    make_dot_gradf(grad::Function, x = nothing; w = similar(x), Tf::DataType)

Return a function with arguments `(v,x)`
for computing a dot product between an array `v`
and the gradient `grad(x)` of some real-valued function `f(x)`,
evaluated at array `x`,
typically for use in a line-search method.

- For a single-argument gradient function `grad(x)`,
  this returns the version
  `(v, x) -> dot(v, grad(x))`.
  This version will be allocating,
  unless `grad` has its own internal workspace due to a closure.
  This version ignores all the arguments except `grad`.

- For a two-argument mutating gradient function `grad!(w, x)`,
  this returns a function
  `(v, x) -> dot(v, grad!(w, x))`
   using the keyword argument `w` as the work array.

If `f(x)` maps an Array `x`
of elements with units `Ux`
into real numbers with units `Uf`,
then the gradient `∇f` has units `Uf/Ux`.
Those units are relevant to defining the work array `w`.

# in
- `grad::Function` see above
- `x` an array whose `size` and `eltype` is used to allocate `w`
  (or `nothing`)

# option
- `Tf::DataType = typeof(one(eltype(x)))`
  Specify `eltype` of function `f(x)`, defaulting to unitless.
- `w = similar(x,  typeof(oneunit(Tf) / oneunit(eltype(x))))`
  work space for gradient calculation, with appropriate units (if needed).

# out
- `(v, x) -> dot(v, grad([w,] x))` or equivalent.

Example.
Consider the cost function
`f(x) = sum(sin.(x))`
that has gradient
`g(x) = cos.(x)`.
The simple way here is
`make_dot_gradf(g)`,
which will return the function
`(v, x) = dot(v, g(x)) === dot(v, cos.(x))`.
That will work fine,
but the `cos.(x)` step will be allocating.
A better way (for repeated use) is
```julia
function g!(work, x)
   @. work = cos(x) # mutating
   return work
end
```
Then
`make_dot_gradf(g!, x)`
will return a function `(v,x)`
that does not allocate
(except when first generated
via a closure).
For this particular example,
an even better approach
is to directly define
`dot_gradf(v,x) = sum(vx -> vx[1]*cos(vx[2]), zip(v,x))`,
but `make_dot_gradf` is here
for more general cases.
"""
function make_dot_gradf(
    g!::Function,
    x::AbstractArray{Tx},
    ;
    Tf::DataType = typeof(one(Tx)),
    w = similar(x, typeof(oneunit(Tf) / oneunit(Tx))),
) where {Tx <: Number}
    _narg(g!) == 2 || error("g! needs '(w,x)' arguments for mutating version!")
    return (v, x) -> dot(v, g!(w, x))
end

function make_dot_gradf(∇f::Function)
    narg = _narg(∇f)
    narg == 2 && error("need 'x' argument for mutating version!")
    narg == 1 || error("∇f needs one argument")
    return (v, x) -> dot(v, ∇f(x))
end

#=
# type inference issues:
function make_dot_gradf(
    grad::Function,
    x = nothing,
    ;
    Tf::DataType = isnothing(x) ? Nothing : typeof(one(eltype(x))),
    g = isnothing(x) ? nothing : similar(x, typeof(oneunit(Tf) / oneunit(eltype(x)))),
)
    narg = first(methods(grad)).nargs - 1
    narg ∈ (1,2) || error("bad narg=$narg")
    narg == 1 && (!isnothing(x) || !isnothing(g)) && @warn("ignoring x,g")
    narg == 2 && (isnothing(x) || isnothing(g) || axes(x) != axes(g) || eltype(x) != eltype(g)) &&
        error("grad!(w,f) requires 'x' argument")
    return narg == 1 ? (v, x) -> dot(v, grad(x)) : (v, x) -> dot(v, grad(w, x))
end
=#


# handle the case where cf returns a constant
make_dot_curvf(curv::RealU) = (v, x) -> curv * sum(abs2, v)


"""
    make_dot_curvf(curv::Function, x; c = similar(x))

Make a function `(v,x)`
that computes the dot product between `abs2.(v)` and `curv(x)`
where `curv(x)` is a curvature function
associated with a quadratic majorizer
for some cost function `f(x)`,
typically for use in a line-search method.

- For a single-argument gradient function `curv(x)`,
  this returns the version equivalent to
  `(v, x) -> dot(abs2.(v), curv(x))`.
  This version will be allocating,
  unless `curv` has its own internal workspace due to a closure.
  This version ignores all the arguments except `curv`.

- For a two-argument mutating gradient function `curv!(c, x)`,
  this returns a function equivalent to
  `(v, x) -> dot(abs2.(v), curv!(c, x))`
   using the keyword argument `c` as the work array.

- If `curv` is simply a real number (a Lipschitz constant),
  possibly with units,
  then this returns `(v, x) -> curv * sum(abs2, v)`.

If `f(x)` maps an Array `x`
of elements with units `Ux`
into real numbers with units `Uf`,
then its curvature has units `Uf/Ux^2`.
Those units are relevant to defining the work array `w`.

# in
- `curv::Function` see above
- `x` an array whose `size` and `eltype` is used to allocate `w`
  (or `nothing`)

# option
- `Tf::DataType = typeof(one(eltype(x)))`
  Specify `eltype` of function `f(x)`, defaulting to unitless.
- `w = similar(x,  typeof(oneunit(Tf) / oneunit(eltype(x))^2))`
  work space for gradient calculation, with appropriate units (if needed).

# out
- `(v, x) -> dot(abs2.(v), curv([w,] x))` or equivalent.
"""
function make_dot_curvf(
    cf!::Function,
    x::AbstractArray{Tx},
    ;
    Tf::DataType = typeof(one(Tx)),
    w = similar(x, typeof(oneunit(Tf) / oneunit(Tx)^2)),
) where {Tx <: Number}
    _narg(cf!) == 2 || error("cf! needs '(w,x)' arguments for mutating version!")
#   return (v, x) -> dot(v, cf!(w, x))
    return (v, x) -> sum(vc -> abs2(vc[1]) * vc[2], zip(v, cf!(w, x)))
end

function make_dot_curvf(cf::Function)
    narg = _narg(cf)
    narg == 2 && error("need 'x' argument for mutating version!")
    narg == 1 || error("cf needs one argument")
#   return (v, x) -> dot(abs2.(v), cf(x))
    return (v, x) -> sum(vc -> abs2(vc[1]) * vc[2], zip(v, cf(x)))
end


# this will fail (as it should) if the units of uj and vj are incompatible
_ls_mm_worktype(uj, vj, α::Real) =
    typeof(oneunit(eltype(uj)) + α * oneunit(eltype(vj)))
#=
    promote_type(
        eltype(u[j]),
        eltype(v[j]),
        typeof(oneunit(eltype(u[j])) + α * oneunit(eltype(u[j]))),
    )
=#


"""
    LineSearchMMWork{Tz <: AbstractVector{<:AbstractArray}}

Workspace for storing ``z_j = u_j + α v_j`` in MM-based line search.

If all of those ``z_j`` arrays had the same `eltype`,
then we could save memory
by allocating just the longest vector needed.
But for Unitful data they could have different eltypes and sizes,
which would require a lot of ``reinterpret`` and ``reshape`` to handle.
So for now it is easier
to just allocate separate work arrays for each `j`.
"""
mutable struct LineSearchMMWork{Tz <: AbstractVector{<:AbstractArray}}
    zz::AbstractVector{<:AbstractArray}

    function LineSearchMMWork(
        uu::AbstractVector{<:AbstractArray},
        vv::AbstractVector{<:AbstractArray},
        α::Real,
    )

        axes(uu) == axes(vv) || error("incompatible u,v axes")
        all(j -> axes(uu[j]) == axes(vv[j]), eachindex(uu)) ||
            error("incompatible uj,vj axes")

        zz = [Array{_ls_mm_worktype(uu[j], vv[j], α)}(undef, size(uu[j])) for j in eachindex(uu)]
        Tz = typeof(zz)
        return new{Tz}(zz)
    end
end


"""
    LineSearchMMState{...}

Mutable struct for MM-based line searches.
"""
mutable struct LineSearchMMState{
    Tu <: AbstractVector{<:AbstractArray},
    Tv <: AbstractVector{<:AbstractArray},
    Tg <: AbstractVector{<:Function},
    Tc <: AbstractVector{<:Union{Function,RealU}},
#   Tw <: AbstractVector{<:AbstractArray},
    Tw <: LineSearchMMWork,
    Tα <: Real,
}
    uu::Tu # vector of B_j x	
    vv::Tv # vector of B_j d	
    dot_gradf::Tg # vector of <z,∇f> functions
    dot_curvf::Tc # vector of <|z|²,ω_f> functions
    work::Tw
    α::Tα

    function LineSearchMMState(
        uu::Tu,
        vv::Tv,
        dot_gradf::Tg,
        dot_curvf::Tc,
        α::Tα = 0f0,
        work::Tw = LineSearchMMWork(uu, vv, α),
    ) where {
        Tu <: AbstractVector{<:AbstractArray},
        Tv <: AbstractVector{<:AbstractArray},
        Tg <: AbstractVector{<:Function},
        Tc <: AbstractVector{<:Union{Function,RealU}},
        Tw <: LineSearchMMWork,
        Tα <: Real,
    }

        all(==(axes(work.zz)), axes.((uu, vv, dot_gradf, dot_curvf))) ||
            error("incompatible axes")
        all(j -> axes(uu[j]) == axes(vv[j]), eachindex(uu)) ||
            error("incompatible u,v axes")
        all(j -> _ls_mm_worktype(uu[j], vv[j], α) == eltype(work.zz[j]), eachindex(uu)) ||
            error("incompatible work type")

        new{Tu, Tv, Tg, Tc, Tw, Tα}(uu, vv, dot_gradf, dot_curvf, work, α)
    end
end


# todo: iterate

# MM-based line search for step size α
# using h(α) = sum_j f_j(uj + α vj)
# \dot{h}(α) = sum_j v_j' * ∇f_j(u_j + α v_j)
function _update!(state::LineSearchMMState)
    uu = state.uu
    vv = state.vv
    zz = state.work.zz
    α = state.α
    dot_gradf = state.dot_gradf
    dot_curvf = state.dot_curvf

    for j in eachindex(zz) # todo: parfor
        @. zz[j] = uu[j] + α * vv[j]
    end
    derh = sum(j -> real(dot_gradf[j](vv[j], zz[j])), eachindex(zz))
    curv = sum(j -> dot_curvf[j](vv[j], zz[j]), eachindex(zz))
    curv < 0 && error("bug: curv=$curv < 0")
    if curv > 0
        state.α -= derh / curv
    end
end


"""
    α = line_search_mm(u, v, gradf, curvf; α0 ...)

Line-search based on majorize-minimize (MM) approach
for a very general family of 1D cost functions of the form
``h(α) = \\sum_{j=1}^J f_j(u_j + α v_j)``
where each function ``f_j(t)`` has a quadratic majorizer of the form
```math
q_j(t;s) = f_j(t) + \\nabla f_j(s) (t - s) + 1/2 \\|t - s\\|^2_{C_j(s)}
```
where ``C_j(⋅)`` is diagonal matrix of curvatures.
(It suffices for each ``f_j`` to have a Lipschitz smooth gradient.)

Each function ``f_j : \\mathcal{X}_j ↦ \\mathbb{R}``
where conceptually
``\\mathcal{X}_j ⊆ \\mathbb{R}^{M_j}``,
but we allow more general domains.

# in
- `u` vector of ``J`` arrays ``u_1,…,u_J`` (typically vectors)
- `v` vector of ``J`` arrays ``v_1,…,v_J`` (typically vectors)
  We require `axes(u_j) == axes(v_j)` for all ``j=1,…J``.
- `gradf` vector of ``J`` functions return gradients of ``f_1,…,f_J``
- `curvf` vector of ``J`` functions `z -> curv(z)` that return a scalar
  or a vector of curvature values for each element of ``z``

# option
- `α0::Real = 0f0` initial guess
- `ninner::Int = 5` # max number of inner iterations of MM line search
- `fun` User-defined function to be evaluated with the `state`
  initially and then after each iteration.
- `out::Union{Nothing,Vector{Any}} = nothing`
  optional place to store result of `fun`:
  `fun(state,0), fun(state,1), ..., fun(state,ninner))`.
   (All `missing by default.) This is a vector of length `ninner+1`.
- `work = LineSearchMMWork(u, v, α)` pre-allocated work space for ``u_j+α v_j``
- `dot_gradf::AbstractVector{<:Function} = make_dot_gradf.(gradf)`
  See `make_dot_gradf`.
- `dot_curvf::AbstractVector{<:Function} = make_dot_curvf.(curvf)`
  See `make_dot_curvf`.

# output
- `α` final iterate

This function mutates the optional arguments `out` and `work`.
"""
function line_search_mm(
    uu::AbstractVector{<:Any},
    vv::AbstractVector{<:Any},
    gradf::AbstractVector{<:Function}, # ignored if dot_gradf provided
    curvf::AbstractVector{<:Union{Function,Number}}, # ignored if dot_curvf provided
    ;
    α0::Real = 0f0,
    work::Tw = LineSearchMMWork(uu, vv, α0),
    out::Union{Nothing,Vector{Any}} = nothing,
    ninner::Int = 5,
    fun::Function = (state,iter) -> missing,
    dot_gradf::AbstractVector{<:Function} = make_dot_gradf.(gradf),
    dot_curvf::AbstractVector{<:Function} = make_dot_curvf.(curvf),
) where {Tw <: LineSearchMMWork}

    !isnothing(out) && length(out) < ninner+1 && throw("length(out) < $(ninner+1)")

    state = LineSearchMMState(uu, vv, dot_gradf, dot_curvf, α0, work)
    if !isnothing(out)
        out[1] = fun(state, 0)
    end

    for iter in 1:ninner
        _update!(state)
        if !isnothing(out)
             out[iter+1] = fun(state, iter)
        end
        state.α == 0 && break
    end

    return state.α
end


#=
todo: special case of a single function

"""
    (x,out) = ncg(grad, curv, x0, ...)

Special case of `ncg` (nonlinear CG) for minimizing a cost function
whose gradient is `grad(x)`
and that has a quadratic majorizer with diagonal Hessian given by
`curv(x)`.
Typically `curv = (x) -> L` where `L` is the Lipschitz constant of `grad`.
"""
function ncg(
    grad::Function,
    curv::Function,
    x0::AbstractArray{<:Number} ;
    kwargs...,
)

    return ncg([I], [grad], [curv], x0; kwargs...)
end
=#
