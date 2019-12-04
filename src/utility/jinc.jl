#=
jinc.jl
=#

export jinc

using SpecialFunctions: besselj1 #grab bessel function
using Plots: plot
using Test: @inferred

"""
    jinc(x)
return `jinc(x) = J1(pi*x)/(2x)`, where `J1` is a Bessel function of the first kind.
units of `x` are typically cycles/m
"""
function jinc(x::Real)
    if (x == 0) return pi/4 end
    y = abs(x)
    besselj1(pi*y) / (2*y)
end

"""
    jinc(:test)
self test
"""
function jinc(x::Symbol)
    x != :test && throw("non-test, symbolic input to jinc")
    r = LinRange(-10,10,201)
    y = @inferred jinc.(r)
    plot(r, y)
    true
end
