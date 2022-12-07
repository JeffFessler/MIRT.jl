# test/regularize/diffs.jl

using MIRT: diff_map, diff_forw, diff_adj
import MIRT: diff2d_forw, diff2d_adj, diff_length, diff2d_map # deprecated

using LinearMapsAA: LinearMapAA
using Test: @test, @test_throws


@test_throws ArgumentError diff_forw(ones(3), dims=2)
@test_throws ArgumentError diff_forw(ones(3), dims=[1,2])
@test_throws ArgumentError diff_forw(ones(1,2,1), dims=(1,2))
@test diff_forw(ones(2,4,6) ; dims=2) == zeros(2*3*6)

for N in [(3,), (3,4), (2,3,4), (4,3,2,2)]
    T = diff_map(N)
    @test Matrix(T)' == Matrix(T')
    @test T.name == "diff_map"
    T = diff_map(N, dims=[1])
    @test Matrix(T)' == Matrix(T')
    if length(N) >= 2
        for dims in [2, (2,), [1,2], (1,2)]
            T = diff_map(N, dims=dims)
            @test Matrix(T)' == Matrix(T')
        end
    end
    if length(N) >= 3
        for dims in [3, (1,3), [2,3], (1,2,3)]
            T = diff_map(N, dims=dims)
            @test Matrix(T)' == Matrix(T')
        end
    end
end

if true # test old 2D versions
    N = (3,5)
    x = rand(N...)
    @test diff2d_forw(x) == diff_forw(x)
    d = rand(sum(diff_length(N,dim) for dim=1:2))
    @test diff2d_adj(d, N... ; out2d=false) == vec(diff_adj(d, N))
    @test diff2d_adj(d, N... ; out2d=true) == diff_adj(d, N)
    @test Matrix(diff2d_map(N...)) == Matrix(diff_map(N))
end

N = (1,2)
@test_throws ArgumentError T = diff_map(N)
@test_throws ArgumentError T = diff_map(N, dims=1)
T = diff_map(N, dims=2)
@test Matrix(T)' == Matrix(T')
