#=
downsample.jl
timing tests
=#

import MIRT
using BenchmarkTools: @btime


function downsample3_time( ; T=Float32, dims=(100,200,40), down=(2,2,2))
    x = ones(T, dims...)
    fp(x) = MIRT.downsample3_perm(x, down)
    fl(x) = MIRT.downsample3_loop(x, down; T)
    @assert fp(x) == fl(x)
    @btime $fp($x) # fastest
    @btime $fl($x) # 2x slower, many more allocations
    nothing
end

#= ir28 with 1.4.2
  3.3 ms (63 allocations: 8.01 MiB)
  5.4 ms (200011 allocations: 9.54 MiB)
=#

downsample3_time()
