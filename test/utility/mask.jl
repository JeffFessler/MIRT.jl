# test/utility/mask.jl

using MIRT: embed, embed!, mask_or, mask_outline, maskit, getindex!

using Test: @test, @testset, @test_throws, @inferred
#using SparseArrays: sparse


function embed_test()
    mask = [false true true; true false false]
    a = [0 2 3; 1 0 0]; v = a[mask]
    b = similar(a)
    @test (@inferred embed(v,mask)) == a
    @test (@inferred embed([v 2v], mask)) == cat(dims=3, a, 2a)
    @test embed!(b, v, mask ; filler=-1) == a + (-1) * .!mask
#   @test embed(sparse(1:3),mask) == sparse([0 2 3; 1 0 0]) # later
    true
end


function getindex!_test()
    N = (2^2,2^3)
    mask = rand(N...) .< 0.5
    K = sum(mask)
    T = ComplexF32
    x = randn(T, N...)
    y0 = x[mask]
    y1 = Array{T}(undef, K)
    @test 0 == @allocated getindex!(y1, x, mask)
    @test y0 == y1
    true
end


function mask_or_test()
    mask2 = trues(3,4)
    @test (@inferred mask_or(mask2)) === mask2
    @test sum(mask_or(mask2)) == length(mask2)
    mask3 = trues(3,4,5)
    @test (@inferred mask_or(mask3)) == trues(3,4)
    @test_throws String mask_or(trues(1,))
    true
end


function mask_outline_test()
    mask2 = trues(3,4)
    @test (@inferred mask_outline(mask2)) == falses(3,4)
    mask3 = trues(3,4,5)
    @test (@inferred mask_outline(mask3)) == falses(3,4)
    true
end


function maskit_test()
    mask = [false true true; true false false]
    @inferred maskit([7 2 3; 1 7 7], mask)
    data = [7 2 3; 1 7 7]
    @test maskit(data, mask) == [1,2,3]
    data = repeat(data, 1, 1, 2)
    @test maskit(data, mask) == repeat(1:3, 1, 2)
    @test_throws DimensionMismatch maskit([0], mask)
    true
end


@testset "embed" begin
    @test embed_test()
end
@testset "getindex!" begin
    @test getindex!_test()
end
@testset "mask_or" begin
    @test mask_or_test()
end
@testset "mask_outline" begin
    @test mask_outline_test()
end
@testset "maskit" begin
    @test maskit_test()
end
