# test/runtests.jl

using Test
using MIRT

include("../../data/z-list.jl")

include("../algorithm/z-list.jl")
include("../fbp/z-list.jl")
include("../io/z-list.jl")
include("../plot/z-list.jl")
include("../regularize/z-list.jl")
include("../system/z-list.jl")
