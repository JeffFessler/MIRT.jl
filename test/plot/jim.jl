# jim.jl

using MIRT: jim
using LaTeXStrings
using Test: @test, @test_throws


jim()
jim(:keys)
jim(:clim)
@test typeof(jim(:defs)) <: Dict

@test_throws String jim(:bad)

jim(ones(4,3), title="test2", xlabel=L"x")
jim(rand(4,3,5), title=L"test3 x^2_i")
jim(1:4, 5:9, zeros(4,5), title="test3", ylabel=L"y")
jim(zeros(4,5), x=1:4, y=5:9, title="test3")
jim(rand(6,4), fft0=true)
jim(x=1:4, y=5:9, rand(4,5), title="test4")
jim(x=-9:9, y=9:-1:-9, (-9:9) * (abs.((9:-1:-9) .- 5) .< 3)', title="rev")
jim(ones(3,3)) # uniform
jim(:abswarn, false)
jim(complex(rand(4,3)))
jim(complex(rand(4,3)), "complex")
jim(:abswarn, true)
jim(rand(4,5), color=:hsv)
