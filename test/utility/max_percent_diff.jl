# test/utility/max_percent_diff.jl

using MIRT: max_percent_diff

using Test: @test


@test max_percent_diff([200], [202]) ≈ 100/101
@test max_percent_diff([0 200], [0 202]) ≈ 100/101
@test max_percent_diff([0 200], [0 198], maxboth=true) ≈ 1
#≈ does not work with 0 values (it's proportional to the second arg)
@test max_percent_diff([0 100], [0 200], normalize=true) < .01
