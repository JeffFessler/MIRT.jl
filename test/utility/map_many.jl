# map_many.jl

using MIRT: map_many

using Test: @test, @inferred

N = (4,3)
fun = x -> cat(dims=3, x, 4x)
x1 = ones(N)
x = cat(dims=3, x1, 2x1)

#@inferred map_many(fun, x, N) # fails because fun could be anything!
@test map_many(fun, x, N) == cat(dims=4, fun(x1), fun(2x1))
