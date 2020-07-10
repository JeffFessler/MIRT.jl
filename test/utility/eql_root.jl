# eql_root.jl

using MIRT: eql_root

using Test: @test, @test_throws, @inferred


tests = [
1 1 -1
0 1 2
4 5 0
1 0 9 # 3 both times
# 1 0 -9; # -3 both times. omitted because out of prog scope
4 4 5
] # (2x+5)(2x-1), so 1/2.

predicted = [-1, 1, 0, 3, 1/2]
results = eql_root(tests[:,1], tests[:,2], tests[:,3])
  
@test results == predicted
@test eql_root(1, 1, -1) == -1
