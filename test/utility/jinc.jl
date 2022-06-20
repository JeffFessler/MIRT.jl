# jinc.jl

using MIRT: jinc

#using Plots: plot
using Test: @test, @inferred

@test jinc(0) ≈ pi/4

jinx = x -> jinc.(x)
r = LinRange(-10,10,201)
y = @inferred jinx(r)

for T in (Int, Float16, Float32, Float64)
    @inferred jinc(T(0))
    @inferred jinc(T(1))
end

#plot(r, y)
