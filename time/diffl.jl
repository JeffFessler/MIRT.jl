#=
diffl.jl
timing tests
=#

using MIRT: diffl!, diff_forw
using Test: @test
using BenchmarkTools: @btime


#=
# timing tests show ~40% speedup for in-place version
=#


# time vs diff
function time0()
	N = (2^10, 2^9)
    x = rand(N...)
    g0 = rand(N...)
    g1 = rand(N...)
	g0 = diff(x, dims=2)
    diffl!(g1, x, 2)
	@test g1[:,2:end] == g0
    @btime diffl!($g1, $x, 2)		# 187 528 b
    @btime g0 = diff($x, dims=2)	# 305 3.99 MiB
    @btime g1 = diffl($x, 2)		# 305 4.00 MiB
	nothing
end

function time1()
    N = 2^10
    x = rand(N)
    g = rand(N)
    @btime diffl!($g, $x)
    @btime diffl!($g, $x ; add=true)
    @btime g = diffl($x)
    nothing
end

function time2()
	N = (2^10, 2^9)
    x = rand(N...)
    g = rand(N...)
    @btime diffl!($g, $x, 1)
    @btime g = diffl($x)
    nothing
end

using MIRT: diff_forw, diff_adj

# in-place version diffl! about 40% faster than diff_forw
function time3()
	N = (2^5, 2^6, 2^7)
    x = rand(N...)
	d = 2
    g2 = diffl(x, d)
    h2 = diff_forw(x, dims=d)
	@test vec(g2[:,2:end,:]) == h2
    diffl!(g2, x, d)
	@test vec(g2[:,2:end,:]) == h2
    @btime diffl!($g2, $x, $d)
    @btime diff_forw($x, dims=$d)
end

#time0()
