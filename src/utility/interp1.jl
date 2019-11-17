#=
interp1.jl
emulate matlab's interp1()
=#

export interp1

using Plots: scatter, plot!, gui
using Interpolations
using Test: @test, @test_throws, @inferred


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
		extrap::Any = 0)
#		extrap::Union{<:Interpolations.BoundaryCondition, <:Number} = 0 # fails?

	fun = interpolate((x,), y, how)
	if !isnothing(extrap)
		fun = extrapolate(fun, extrap)
	end
	fun.(xi)
end


"""
    interp1(:test)
self test
"""
function interp1(test::Symbol)
	@test (@inferred interp1(1:4, 2:5, 1.5)) == 2.5
	@test (@inferred interp1(1:2, [3im,4im], 1.5)) == 3.5im
	x = LinRange(-1, 1, 21)
	y = sin.(x*π)
	xi = LinRange(-1, 1, 101) * 2
	@test_throws BoundsError interp1(x, y, xi, extrap=nothing)
	@test_throws BoundsError interp1(x, y, xi, extrap=Throw())
	scatter(x, y, label="")
	y0 = @inferred interp1(x, y, xi) # 0 extrapolation
	plot!(xi, y0, label="0", ylim=[-2,2])
	for ex in (Flat(), Line(), Periodic(), Reflect())
		yy = @inferred interp1(x, y, xi, extrap=ex)
		plot!(xi, yy, label="$ex")
	end
	gui()
	true
end
