# eql_root.jl

export eql_root


"""
    x = eql_root(a,b,c)

Numerically stable method for computing the positive root
of the quadratic polynomial `-ax^2 - 2bx + c, a >= 0`.
Assumes solvable equations; will throw otherwise.

# in
- `a` The negative of the `x^2` term. Must be positive.
- `b` Half the negative of the `x` term.
- `c` The constant term.

# out
- `x` The positive root that satisfies `0 = -ax^2 - 2bx + c`.
"""
function eql_root(a::Real, b::Real, c::Real)
    (a < 0) && throw(DomainError(a, "a must be enonnegative"))

    if a == 0
        return c / b / 2 # trivially solve case where a == 0
    end

    det = sqrt(b^2 + a * c) # determinant / 2

    # handle positive and negative b appropriately:
    return (b > 0) ? c / (det + b) : (det - b) / a
end


function eql_root(a, b, c)
    (size(a) != size(b) || size(b) != size(c)) && throw(DimensionMismatch("all arguments must share dimensions"))

    eql_root.(a, b, c)
end
