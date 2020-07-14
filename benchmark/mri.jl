# exp_xform
SUITE["exp_xform"] = BenchmarkGroup()
let T=ComplexF32, N=100, M=10000, D=20, L=10
    rng = MersenneTwister(0)
    X = randn(rng, T, N, L)
    U = randn(rng, T, D, N)
    V = randn(rng, T, D, M)

    for mode in (:matrix, :row, :element, :column)
        # timing tests: :matrix is fastest, with :row a close 2nd
        fn(x, u, v) = exp_xform(x, u, v; mode=mode)
        SUITE["exp_xform"][string(mode)] = @benchmarkable $fn($X, $U, $V)
    end
end
