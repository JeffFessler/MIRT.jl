#=
interp1.jl
emulate matlab's interp1()
=#

export interp1

using Interpolations


"""
    yi = interp1(x, y, xi ; how=Gridded(Linear()), extrap=0)

1D interpolation of `y = f(x)` at points `xi`

In:
- `x::AbstractVector{<:Real}`
- `y::AbstractVector{<:Number}`

Option:
- `how::Interpolations.InterpolationType` default `Gridded(Linear())`
- `extrap::Any` how to extrapolate, e.g., `Flat()`; default `0`
other options from Interpolations.jl are `Line()` `Periodic()` `Reflect()` `Throw()`

Output is same size as input `xi`
"""
function interp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Number}, xi ;
    how::Interpolations.InterpolationType = Gridded(Linear()),
    extrap::Any = 0,
#   extrap::Union{<:Interpolations.BoundaryCondition, <:Number} = 0 # fails?
)

    fun = interpolate((x,), y, how)
    if !isnothing(extrap)
        fun = extrapolate(fun, extrap)
    end
    fun.(xi)
end
