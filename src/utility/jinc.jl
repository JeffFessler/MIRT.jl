#=
jinc.jl
=#

export jinc

using SpecialFunctions: besselj1 # grab bessel function


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
