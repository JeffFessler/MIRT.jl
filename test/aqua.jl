using MIRT: MIRT
import Aqua
using Test: @testset

@testset "aqua" begin
    Aqua.test_all(MIRT; stale_deps = (ignore = [:AVSfldIO],))
end
