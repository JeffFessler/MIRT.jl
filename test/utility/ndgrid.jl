# ndgrid.jl

using MIRT: ndgrid

using Test: @test, @inferred


@inferred ndgrid(1:3, 2:4)
@test ndgrid(1:3, 2:4) == ([1 1 1; 2 2 2; 3 3 3], [2 3 4; 2 3 4; 2 3 4])
@test size(ndgrid(1:3, 4:5, 6:9)[3]) == (3,2,4)
