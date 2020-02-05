# eql_root.jl

export eql_root

using Test: @test, @test_throws, @inferred
using BenchmarkTools: @btime

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
function eql_root(a, b, c)
    any(a .< 0) && throw(DomainError(a,"a must be entirely nonnegative"))
    (size(a) != size(b) || size(b) != size(c)) && throw(DimensionMismatch("all arguments must share dimensions"))
    
    T = promote_type(eltype(a), eltype(b), eltype(c), Float32)
    x = zeros(T,size(a))
    
    j = (a .== 0)
    x[j] .= c[j] ./ b[j] / 2 # trivially solve equations where a == 0
    det = sqrt.(b.^2 + a .* c) # determinant / 2
    # check for equations where a is positive and b is positive
    j = (a .> 0) .& (b .> 0)
    x[j] .= c[j] ./ (det[j] + b[j])
    # check for equations where a is pos and b is nonpos
    j = (a .> 0) .& (b .<= 0)
    x[j] .= (det[j] - b[j]) ./ a[j]

    return x
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
    return true
end
