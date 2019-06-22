#= mask.jl
Methods related to a image support mask:
	embed, maskit, mask_or, mask_outline
2019-06-22 Jeff Fessler, University of Michigan
=#

using Test
using ImageFiltering: imfilter, centered
#using DSP: conv2
#using SparseArrays: sparse, findnz, AbstractSparseVector


"""
`mask_or(mask)`
compress 3D mask to 2D by logical `or` along `z` direction
"""
function mask_or(mask::AbstractArray{Bool})
	ndim = ndims(mask)
	return	ndim == 2 ? mask :
			ndim == 3 ? dropdims(sum(mask, dims=3) .> 0, dims=3) :
			throw("ndims(mask) = $ndim")
end


"""
`mask_or(:test)`
self test
"""
function mask_or(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask2 = trues(3,4)
	@test mask_or(mask2) == mask2
	mask3 = trues(3,4,5)
	@test mask_or(mask3) == trues(3,4)
	true
end


"""
`mask_outline(mask)`
return outer boundary of 2D mask (or mask_or for 3D)
"""
function mask_outline(mask::AbstractArray{Bool})
	mask2 = mask_or(mask)
	tmp = imfilter(mask2, centered(ones(Int32,3,3))) # dilate
#	tmp = conv2(Float32.(mask2), ones(Float32,3,3))
#	tmp = tmp[2:end-1,2:end-1] # 'same'
	return (tmp .> 1) .& (.! mask2)
end


"""
`mask_outline(:test)`
self test
"""
function mask_outline(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask2 = trues(3,4)
	@test mask_outline(mask2) == falses(3,4)
	true
end


"""
`array = embed(v, mask)`

embed vector `v` of length `sum(mask)`
into elements of an array where `mask` is `true`
"""
function embed(v::AbstractVector{<:Number}, mask::AbstractArray{Bool})
	array = zeros(eltype(v), size(mask))
	array[mask] .= v
	return array
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
`embed(:test)`
"""
function embed(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask = [false true true; true false false]
	@test embed(1:3,mask) == [0 2 3; 1 0 0]
#	@test embed(sparse(1:3),mask) == sparse([0 2 3; 1 0 0]) # later
	true
end


"""
`maskit(x::AbstractArray{<:Number})`
opposite of embed
"""
function maskit(x::AbstractArray{<:Number}, mask::Array{Bool})
	dim = size(x)
	x = reshape(x, length(mask), :)
	x = x[mask[:],:] # reshape(mask, prod(_dim()))]
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
`maskit(:test)`
"""
function maskit(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	mask = [false true true; true false false]
	@test maskit([7 2 3; 1 7 7], mask) == [1,2,3]
	true
end


#=
@test mask_or(:test)
@test mask_outline(:test)
@test embed(:test)
@test maskit(:test)
=#
