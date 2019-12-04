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
return type is `promote_type(typeof(x), Float32)`
"""
function jinc(x::Real)
    T = promote_type(typeof(x), Float32)
    if (x == 0) return convert(T, pi/4) end
    y = abs(x)
    convert(T, besselj1(pi*y) / (2*y))
end

"""
    jinc(:test) jinc(:plot)
self test, plot
"""
function jinc(test::Symbol)
    if test === :plot
        r = LinRange(-10,10,201)
        return plot(r, jinc.(r))
    end
    !(test === :test) && throw("non-test, symbolic input to jinc")
    jinx = x -> jinc.(x)
    r = LinRange(-10,10,201)
    @inferred jinx(r)
    for T in (Int, Float16, Float32, Float64)
        @inferred jinc(T(0))
        @inferred jinc(T(1))
    end
    jinc(:plot)
    true
end
