#=
exp_xform.jl
timing test
=#

using MIRT: exp_xform
using BenchmarkTools: @btime


# timing test
function exp_xform_time( ; T::DataType = ComplexF32)

    modes = (:element, :row, :column)

    N = 10000
    M = 100
    D = 20
    L = 10
    X = randn(T, N, L)
    U = randn(T, D, N)
    V = randn(T, D, M)

    # timing tests: :matrix is fastest, with :row a close 2nd
    for mode in (:matrix, modes...)
        @info "Case :$mode"
        tmp = (x, u, v) -> exp_xform(x, u, v ; mode=mode)
        @btime $tmp($X, $U, $V)
    end

#=
ir28 with 1.4.2:

:matrix   25.248 ms (8 allocations: 22.90 MiB)
:element 353.762 ms (4000001 allocations: 762.95 MiB)
:row      40.402 ms (1001 allocations: 22.99 MiB)
:column   51.424 ms (70001 allocations: 107.12 MiB)
=#

end

exp_xform_time()
