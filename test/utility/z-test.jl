# test/utility/z-test.jl

using Test: @testset

list = [
"downsample.jl"
"eql_root.jl"
"interp1.jl"
"jinc.jl"
"map_many.jl"
"mask.jl"
"max_percent_diff.jl"
"ndgrid.jl"
"rect.jl"
"reverser.jl"
"rmsd100.jl"
]

for file in list
    @testset "$file" begin
        include(file)
    end
end
