#=
diffl.jl
timing tests
=#

using MIRT: diffl!, diffl_adj!, diffl, diffl_adj # new
using MIRT: diff_forw, diff_adj # old
using Test: @test
using BenchmarkTools: @btime


#=
# timing tests show ~40% speedup for in-place version
=#


function diffl_time()

	@info "time vs diff and vs w/o in-place"

	N = (2^10, 2^9)
    x = rand(N...)
    g0 = rand(N...)
    g1 = rand(N...)
	g0 = diff(x, dims=2)
    diffl!(g1, x, 2)
	@test g1[:,2:end] == g0
    @btime diffl!($g1, $x, 2)		# 210 20 672 b
    @btime g0 = diff($x, dims=2)	# 303  4 3.99 MiB
    @btime g1 = diffl($x, 2)		# 306 22 4.00 MiB


	@info "time vs diff_forw"
	# in-place version diffl! about 40% faster than diff_forw

	N = (2^5, 2^6, 2^7)
    x = rand(N...)
	d = 2
    g2 = diffl(x, d)
    h2 = diff_forw(x, dims=d)
	@test vec(g2[:,2:end,:]) == h2
    diffl!(g2, x, d)
	@test vec(g2[:,2:end,:]) == h2
    @btime diffl!($g2, $x, $d)		# 92 24 800 b
    @btime diff_forw($x, dims=$d)	# 148 6 1.97 MiB


	@info "time adj" # about 40% speed up

	z2 = diff_adj(h2, N ; dims=d)
	z1 = diffl_adj(g2, d)
	@test z1 == z2
	diffl_adj!(z1, g2, d)
	@test z1 == z2

	@btime diffl_adj!(z1, g2, d)		# 361 28 928 b
	@btime diff_adj($h2, $N ; dims=$d)	# 607 70 3.94 MiB

	nothing
end

#diffl_time()
