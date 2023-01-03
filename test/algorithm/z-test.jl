# test/algorithm/z-test.jl

using Test: @testset

list = [
"general/dot-grad.jl"
"general/dot-curv.jl"
"general/ls-mm.jl"
"general/ncg.jl"
"general/ogm_ls.jl"
"general/pogm_restart.jl"
"general/poweriter.jl"
]

for file in list
    @testset "$file" begin
        include(file)
    end
end
