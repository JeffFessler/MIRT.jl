#=
dot-curv.jl
=#

using LinearAlgebra: dot

# type (considering units) for curvature df/dx^2
_curv_type(Tf::Type{<:RealU}, Tx::Type{<:Number}) = 
    typeof(oneunit(Tf) / oneunit(Tx)^2)


"""
    make_dot_curvf(curv::Function, [x, Tf::DataType; w = similar(x)])

Make a function with arguments `(v, x)`
that computes the dot product between
`abs2.(v)` and `curv(x)`
where `curv(x)` is a curvature function
associated with a quadratic majorizer for
some real-valued cost function `f(x)`,
evaluated at array `x`,
typically for use in a line-search method.

- For a single-argument gradient function `curv(x)`,
  this returns the version equivalent to
  `(v, x) -> dot(abs2.(v), curv(x))`.
  This version will be allocating,
  unless `curv`
  has its own internal workspace due to a closure.
  This version expects only the argument
  `curv`.

- For a two-argument mutating gradient function
  `curv!(w, x)`,
  this returns a function
  equivalent to
  `(v, x) -> dot(abs2.(v), curv!(w, x))`
  using the keyword argument `w` as the work array.

- If `curv` is simply a real number (a Lipschitz constant),
  possibly with units,
  then this returns the function
  `(v, x) -> curv * sum(abs2, v)`.

If `f(x)` maps an Array `x`
of elements with units `Ux`
into real numbers with units `Uf`,
then its curvature has units `Uf/Ux^2`.
Those units are relevant to defining the work array `w`.

# in
- `curv::Function` see above
- `x` an array whose `size` and `eltype` is used to allocate `w`

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
    Tf::Type{<:RealU} = typeof(one(real(Tx))),
    ;
    w::AbstractArray{Tw} = similar(x, _curv_type(Tf, Tx)),
) where {Tx <: Number, Tw <: Number}
    _narg(cf!) == 2 || error("cf! needs '(w,x)' arguments for mutating version!")
    axes(x) == axes(w) || error("w,x axes mismatch")
    promote_type(Tw, _curv_type(Tf, Tx)) == Tw ||
        error("w type Tw=$Tw vs $(_curv_type(Tf, Tx))")
    let w = w, cf! = cf!
        (v, x) -> sum(vc -> abs2(vc[1]) * vc[2], zip(v, cf!(w, x)))
    end
end


function make_dot_curvf(cf::Function)
    narg = _narg(cf)
    narg == 2 && error("need 'x' argument for mutating version!")
    narg == 1 || error("cf needs one argument")
#   return (v, x) -> dot(abs2.(v), cf(x))
    return let cf = cf
        (v, x) -> sum(vc -> abs2(vc[1]) * vc[2], zip(v, cf(x)))
    end
end


# case where cf returns a constant
make_dot_curvf(curv::RealU) = let curv = curv
    (v, x) -> curv * sum(abs2, v)
end
