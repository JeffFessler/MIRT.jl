#=
mask.jl
Methods related to a image support mask:
	embed, maskit, mask_or, mask_outline
2019-06-22 Jeff Fessler, University of Michigan
=#

export embed, embed!, mask_or, mask_outline, maskit, mask_test
export getindex!

using Test: @test, @testset, @test_throws, @inferred
using ImageFiltering: imfilter, centered
#using SparseArrays: sparse, findnz, AbstractSparseVector


"""
    getindex!(y::Vector, x::AbstractArray, mask::AbstractArray{Bool})
Equivalent to the in-place `y .= x[mask]`.
"""
@inline function getindex!(y::Vector, x::AbstractArray{T,D},
	mask::AbstractArray{Bool,D},
) where {T,D}
    sum(mask) == length(y) || throw("wrong length")
    count = 1
    @inbounds for i in 1:LinearIndices(mask)[findlast(mask)]
        y[count] = x[i]
        count += mask[i]
    end
    return y
end


"""
    mask_or(mask)
compress 3D mask to 2D by logical `or` along `z` direction
"""
mask_or(mask::AbstractArray{Bool,2}) = mask
mask_or(mask::AbstractArray{Bool,3}) =
	dropdims(sum(mask, dims=3) .> 0, dims=3)
mask_or(mask::AbstractArray{Bool}) = throw("ndims(mask) = $(ndims(mask))")


"""
    mask_or(:test)
self test
"""
function mask_or(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask2 = trues(3,4)
	@test (@inferred mask_or(mask2)) == mask2
	mask3 = trues(3,4,5)
	@test (@inferred mask_or(mask3)) == trues(3,4)
	@test_throws String mask_or(trues(1,))
	true
end


"""
    mask_outline(mask)
return outer boundary of 2D mask (or mask_or for 3D)
"""
function mask_outline(mask::AbstractMatrix{Bool})
	tmp = imfilter(mask, centered(ones(Int32,3,3))) # dilate
#	tmp = tmp[2:end-1,2:end-1] # 'same'
	return (tmp .> 1) .& (.! mask)
end
mask_outline(mask::AbstractArray{Bool,3}) = mask_outline(mask_or(mask))


"""
    mask_outline(:test)
self test
"""
function mask_outline(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask2 = trues(3,4)
	@test (@inferred mask_outline(mask2)) == falses(3,4)
	mask3 = trues(3,4,5)
	@test (@inferred mask_outline(mask3)) == falses(3,4)
	true
end


"""
    embed!(array, v, mask ; filler=0)
embed vector `v` of length `sum(mask)`
into elements of `array` where `mask` is `true`,
setting the remaining elements to `filler` (default 0).
"""
embed!(array::AbstractArray{T,D}, v::AbstractVector{<:Number},
    mask::AbstractArray{Bool,D} ; filler::T = zero(T)
) where {T, D} =
	setindex!(fill!(array, filler), v, mask)


"""
    array = embed(v, mask ; filler=0)

embed vector `v` of length `sum(mask)`
into elements of an array where `mask` is `true`; see `embed!`.
"""
embed(v::AbstractVector{T}, mask::AbstractArray{Bool} ;
	filler::Number = zero(T)) where {T <: Number} =
	embed!(fill(filler, size(mask)), v, mask)


"""
    array = embed(matrix::AbstractMatrix{<:Number}, mask::AbstractArray{Bool})

Embed each column of `matrix` into `mask` then `cat` along next dimension
In:
* `matrix [sum(mask) L]`
* `mask [(N)]`

Out:
* `array [(N) L]`
"""
function embed(x::AbstractMatrix{<:Number}, mask::AbstractArray{Bool})
	L = size(x,2)
	out = zeros(eltype(x), prod(size(mask)), L)
	for l=1:L
		out[:,l] = vec(embed(x[:,l], mask))
	end
	reshape(out, size(mask)..., L)
end



#=
this is needed only for "unpacking" sparse system matrices so ignore for now
# image_geom_embed_sparse
# called by image_geom_embed
# function _embed_sparse(x::Array{T}) where {T <: Number}
function embed(x::AbstractSparseVector{<:Number},
				mask::AbstractArray{Bool,N} where N)
	i, v = findnz(x)
	ind = findall(mask)
	j = ind(j)
	return sparsevec(i, j, a, size(x,1), length(mask))
end
=#


"""
    embed(:test)
"""
function embed(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask = [false true true; true false false]
	a = [0 2 3; 1 0 0]; v = a[mask]
	b = similar(a)
	@test (@inferred embed(v,mask)) == a
	@test (@inferred embed([v 2v], mask)) == cat(dims=3, a, 2a)
	@test embed!(b, v, mask ; filler=-1) == a + (-1) * .!mask
#	@test embed(sparse(1:3),mask) == sparse([0 2 3; 1 0 0]) # later
	true
end


"""
    maskit(x::AbstractArray{<:Number})
opposite of embed
"""
function maskit(x::AbstractArray{<:Number}, mask::Array{Bool})
	dim = size(x)
	x = reshape(x, length(mask), :)
	x = x[vec(mask),:] # reshape(mask, prod(_dim()))]
	if length(dim) == ndims(mask)
		x = dropdims(x, dims=2) # squeeze
	elseif length(dim) > ndims(mask)
		x = reshape(x, :, dim[(1+ndims(mask)):end])
	else
		throw(DimensionMismatch("size(x) = $(size(x))"))
	end
	return x
end


"""
    maskit(:test)
"""
function maskit(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask = [false true true; true false false]
	@inferred maskit([7 2 3; 1 7 7], mask)
	@test maskit([7 2 3; 1 7 7], mask) == [1,2,3]
	true
end


"""
    getindex!(:test)
self test
"""
function getindex!(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	N = (2^2,2^3)
	mask = rand(N...) .< 0.5
	K = sum(mask)
	T = ComplexF32
	x = randn(T, N...)
	y0 = x[mask]
	y1 = Array{T}(undef, K)
	getindex!(y1, x, mask)
	@test y0 == y1
	true
end


"""
    mask_test()
self tests
"""
function mask_test()
	@testset "getindex!" begin
		@test getindex!(:test)
	end
	@testset "mask_or" begin
		@test mask_or(:test)
	end
	@testset "mask_outline" begin
		@test mask_outline(:test)
	end
	@testset "embed" begin
		@test embed(:test)
	end
	@testset "maskit" begin
		@test maskit(:test)
	end
	true
end
