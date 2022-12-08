# dtft.jl

import MIRT
using FFTW: fft
using Random: seed!
using BenchmarkTools: @btime

function dtft_time( ; N::Int=2^10, M::Int=2^11, n_shift::Real=7)
    seed!(0)
    x = randn(ComplexF64, N)
    w = randn(M)

    # time DTFT
    @btime fft($x)
    @btime MIRT.dtft_loop_n($w, $x)
    @btime MIRT.dtft_loop_m($w, $x)
    @btime MIRT.dtft_dist_m($w, $x)
    @btime MIRT.dtft_pmap_m($w, $x)
    @btime MIRT.dtft_matvec($w, $x)

#=
    jf28 results (1.4.2 with 8 threads):

    12.355 μs (49 allocations: 18.97 KiB)
    38.909 ms (5122 allocations: 80.31 MiB)
    42.214 ms (2051 allocations: 32.30 MiB)
    43.330 ms (11949 allocations: 32.64 MiB)
    47.251 ms (18472 allocations: 32.64 MiB)
    56.826 ms (7 allocations: 48.03 MiB)

    conclusion: dtft_loop_n is the fastest
    (also amenable to parfor=pmap but unhelpful!?)
=#
end

dtft_time()
