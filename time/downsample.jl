#=
downsample.jl
timing tests
=#

import MIRT
using BenchmarkTools: @btime

function downsample3_time()
	x = ones(Float32, 100, 200, 40)
	@btime MIRT.downsample3_perm($x, (2, 2, 2)) # fastest
	@btime MIRT.downsample3_loop($x, [2, 2, 2], T=Float32) # 2x slower
	nothing
end

#= ir28 with 1.4.2
  3.3 ms (63 allocations: 8.01 MiB)
  5.4 ms (200011 allocations: 9.54 MiB)
=#

downsample3_time()
