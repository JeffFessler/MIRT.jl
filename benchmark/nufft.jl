# dtft
SUITE["dtft"] = BenchmarkGroup()
let T=ComplexF32, N=2^10, M=2^11, n_shift=7
    rng = MersenneTwister(0)
    x = randn(rng, T, N)
    w = randn(rng, M)
    SUITE["dtft"]["fft"] = @benchmarkable fft($x)
    SUITE["dtft"]["loop_n"] = @benchmarkable MIRT.dtft_loop_n($w, $x)
    SUITE["dtft"]["loop_m"] = @benchmarkable MIRT.dtft_loop_m($w, $x)
    SUITE["dtft"]["dist_m"] = @benchmarkable MIRT.dtft_dist_m($w, $x)
    SUITE["dtft"]["pmap_m"] = @benchmarkable MIRT.dtft_pmap_m($w, $x)
    SUITE["dtft"]["matvec"] = @benchmarkable MIRT.dtft_matvec($w, $x)
end
