# test/system/z-test.jl

using Test: @testset

list = [
"Afft.jl"
"Asense.jl"
]

for file in list
    @testset "$file" begin
        include(file)
    end
end
