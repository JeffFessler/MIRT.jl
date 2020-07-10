# interp1.jl

using MIRT: interp1

using Plots: scatter, plot!, gui
using Interpolations: Flat, Line, Periodic, Reflect, Throw
using Test: @test, @test_throws, @inferred


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
