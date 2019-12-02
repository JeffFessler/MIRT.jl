using SpecialFunctions: besselj1 #grab bessel function
using Plots
"""
    jinc(x)
returns jinc(x), defined as J1(x)/x, where J1 is a Bessel function of the first kind.
"""
function jinc(x::Real) #calculates jinc(x)
    if(x == 0) return pi/4 end
    y = abs(x)
    return besselj1(pi*y) / (2*y)
end
function jinc(x::Symbol)
    x != :test && throw("non-test, symbolic input to jinc")
    r = (-10:.01:10)
    plot(r,jinc.(r))
    true
end
