# map_many.jl
# 2019-06-13, Jeff Fessler, University of Michigan

using Test

"""
`y = map_many(fun::Function, x::AbstractArray{<:Any}, idim::Dims)`

apply a function `fun` to leading slices of input `x`;
cousin of `mapslices`

in
* `fun::Function` maps input of size `idim` to output of some size `odim`
* `x` [idim ldim]
out
* `y` [odim ldim]

Example: if `fun` maps array of size (1,2) to array of size (3,4,5)
and if input `x` has size (1,2,7,8)
then output `y` will have size (3,4,5,7,8)
where `y[:,:,:,i,j] = fun(x[:,:,i,j])`
"""
function map_many(fun::Function, x::AbstractArray{<:Any}, idim::Dims)
	xdim = size(x)
	D = length(idim)
	ndims(x) < D && throw(DimensionMismatch("ndims(x) < D=$D)"))
	xdim[1:D] != idim && throw(DimensionMismatch("size(x) vs idim"))

	if xdim == idim
		return fun(x)
	end

	ldim = xdim[(D+1):end]
	L = prod(ldim) # *ldim
	x = reshape(x, prod(idim), :) # [*idim L]
	tmp = fun(reshape((@view x[:,1]), idim)) # [odim]
	odim = size(tmp)
	T = eltype(tmp)
	out = Array{T}(undef, prod(odim), L) # [*odim L]
	out[:,1] = tmp[:]
	for l=2:L
		out[:,l] = fun(reshape((@view x[:,l]), idim))[:]
	end
	return reshape(out, odim..., ldim...) # [odim ldim]
end


"""
`map_many(:test)`
self test
"""
function map_many(test::Symbol)
	N = (4,3)
	fun = x -> cat(dims=3, x, 4x)
	x1 = ones(N)
	x = cat(dims=3, x1, 2x1)

	@test map_many(fun, x, N) == cat(dims=4, fun(x1), fun(2x1))
	true
end


# map_many(:test)
