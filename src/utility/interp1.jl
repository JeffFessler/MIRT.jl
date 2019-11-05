#=
interp1.jl
emulate matlab's interp1()
=#

export interp1

using Plots: scatter, plot!, gui
using Interpolations
using Test: @test, @inferred


"""
    yi = interp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Number}, xi)

1D interpolation of `y = f(x)` at points `xi`

Option:
- `how::Interpolations.InterpolationType` default `Gridded(Linear())`

todo: should have more options for extrapolation etc.

Output is same size as input `xi`
"""
function interp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Number}, xi ;
		how::Interpolations.InterpolationType = Gridded(Linear()))
	fun = interpolate((x,), y, how)
	fun.(xi)
end


"""
    interp1(:test)
self test
"""
function interp1(test::Symbol)
	x = LinRange(0, 2π, 21)
	y = cos.(x)
	xi = LinRange(0, 2π, 101)
	yi = @inferred interp1(x, y, xi)
	@test (@inferred interp1(1:4, 2:5, 1.5)) == 2.5
	@test (@inferred interp1(1:2, [3im,4im], 1.5)) == 3.5im
	scatter(x, y, label="")
	plot!(xi, yi, label="")
	gui()
	true
end

interp1(:test)
