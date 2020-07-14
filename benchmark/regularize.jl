# diff_map
SUITE["diff_map"] = BenchmarkGroup()
let sz=(2^7, 2^7+1)
    x = randn(prod(sz))

    T2d = diff2d_map_old(sz...)
    Tnd = MIRT.diff_map(sz)
    d = T2d * x

    SUITE["diff_map"]["diff_forward"] = BenchmarkGroup()
    SUITE["diff_map"]["diff_forward"]["baseline"] = @benchmarkable $T2d * $x
    SUITE["diff_map"]["diff_forward"]["MIRT"] = @benchmarkable $Tnd * $x

    SUITE["diff_map"]["diff_adjoint"] = BenchmarkGroup()
    SUITE["diff_map"]["diff_adjoint"]["baseline"] = @benchmarkable $T2d' * $d 
    SUITE["diff_map"]["diff_adjoint"]["MIRT"] = @benchmarkable $Tnd' * $d
end

# finite difference
SUITE["finite difference"] = BenchmarkGroup()
let sz = (2^10, 2^9), dims=2
    x = rand(sz...);
    h = diff_forw(x, dims=dims)
    g = diffl(x, dims)
    out = Array{Float64, length(sz)}(undef, sz...)

    SUITE["finite difference"]["diff"] = @benchmarkable diff($x, dims=$dims)
    SUITE["finite difference"]["diffl"] = @benchmarkable diffl($x, $dims)
    SUITE["finite difference"]["diffl!"] = @benchmarkable diffl!($out, $x, $dims)
    SUITE["finite difference"]["diff_forw"] = @benchmarkable diff_forw($x, dims=$dims)

    SUITE["finite difference"]["diff_adj"] = @benchmarkable diff_adj($h, $sz; dims=$dims)
    SUITE["finite difference"]["diffl_adj"] = @benchmarkable diffl_adj($g, $dims)
    SUITE["finite difference"]["diffl_adj!"] = @benchmarkable diffl_adj!($out, $g, $dims)
end
