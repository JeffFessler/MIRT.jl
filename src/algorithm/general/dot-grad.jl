#=
dot-grad.jl
=#

using LinearAlgebra: dot

# type (considering units) for gradient df/dx
_grad_type(Tf::Type{<:RealU}, Tx::Type{<:Number}) =
    typeof(oneunit(Tf) / oneunit(Tx))


"""
    make_dot_gradf(grad::Function, [x, Tf::Type; w = similar(x)])

Make a function with arguments `(v, x)`
that computes the dot product between
an array `v`
and the gradient `grad(x)` of
some real-valued cost function `f(x)`,
evaluated at array `x`,
typically for use in a line-search method.

- For a single-argument gradient function `grad(x)`,
  this returns the version
  `(v, x) -> dot(v, grad(x))`.
  This version will be allocating,
  unless `grad`
  has its own internal workspace due to a closure.
  This version expects only the argument
  `grad`.

- For a two-argument mutating gradient function
  `grad!(w, x)`,
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
- `Tf::Type = typeof(one(eltype(x)))`
  Specify `eltype` of function `f(x)`, defaulting to unitless.

# option
- `w = similar(x, typeof(oneunit(Tf) / oneunit(eltype(x))))`
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
    Tf::Type{<:RealU} = typeof(one(Tx)),
    ;
    w::AbstractArray{Tw} = similar(x, _grad_type(Tf, Tx)),
) where {Tx <: Number, Tw <: Number}
    _narg(g!) == 2 || error("g! needs '(w, x)' arguments for mutating version!")
    axes(x) == axes(w) || error("w,x axes mismatch")
    promote_type(Tw, _grad_type(Tf, Tx)) == Tw || error("w type")
    return (v, x) -> dot(v, g!(w, x))
end


function make_dot_gradf(∇f::Function)
    narg = _narg(∇f)
    narg == 2 && error("need 'x' argument for mutating version!")
    narg == 1 || error("∇f needs one argument")
    return (v, x) -> dot(v, ∇f(x))
end
