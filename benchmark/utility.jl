# downsample
SUITE["downsample"] = BenchmarkGroup()
let T=Float32, sz=(100, 200, 40)
    x = ones(T, sz)
    SUITE["downsample"]["perm"] = @benchmarkable MIRT.downsample3_perm($x, (2, 2, 2))
    SUITE["downsample"]["loop"] = @benchmarkable MIRT.downsample3_loop($x, [2, 2, 2])
end
