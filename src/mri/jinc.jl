using SpecialFunctions #necessary to grab Bessel functions
using Plots
"""
function jinc(x)
returns jinc(x), defined as J1(x)/x, where J1 is a Bessel function of the first kind.
No particular options.
"""
function jinc(x::Real) #calculates jinc(x)
    if(x == 0) return pi/4 end
    y = abs(x)
    return besselj1(pi*y) / (2*y)
end
function jinc(x::Symbol)
    if(x != :test)
        throw("non-test, symbolic input to jinc")
    end
    r = (-10:.01:10)
    plot(r,jinc.(r))
end
