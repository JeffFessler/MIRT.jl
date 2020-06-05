# eql_root.jl

export eql_root

using Test: @test, @test_throws, @inferred


"""
    x = eql_root(a,b,c)
    Numerically stable method for computing the positive root
    of the quadratic polynomial `-ax^2 - 2bx + c, a >= 0`.
    Assumes solvable equations; will throw otherwise.

in
- `a` : The negative of the `x^2` term. Must be positive.
- `b` : Half the negative of the `x` term.
- `c` : The constant term.

out
- `x` : The positive root that satisfies `0 = -ax^2 - 2bx + c`.
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


"""
    eql_root(:test)
self test
"""
function eql_root(x::Symbol)
    x != :test && throw("Must use :test as input for testing")
    tests = [1 1 -1;
             0 1 2;
             4 5 0;
             1 0 9; # 3 both times
             #1 0 -9; # -3 both times. omitted because out of prog scope
             4 4 5] # (2x+5)(2x-1), so 1/2.
    predicted = [-1, 1, 0, 3, 1/2]
    results = eql_root(tests[:,1], tests[:,2], tests[:,3])
  
    @test results == predicted
    @test eql_root(1, 1, -1) == -1
    return true
end
