# test/utility/max_percent_diff.jl

using MIRT: max_percent_diff

using Test: @test, @test_throws

@test max_percent_diff([202], [200]) ≈ 100/101
@test max_percent_diff([0 202], [0 200]) ≈ 100/101
@test max_percent_diff([0 200], [0 198], maxboth=true) ≈ 1
@test max_percent_diff([0 100], [0 200], normalize=true) == 0
@test_throws ErrorException max_percent_diff([0], [2])
@test max_percent_diff([0], [0]) == 0
