# test/io/z-test.jl

using Test: @testset

list = [
"caller_name.jl"
"ir_dump.jl"
"shows.jl"
]

for file in list
    @testset "$file" begin
        include(file)
    end
end
