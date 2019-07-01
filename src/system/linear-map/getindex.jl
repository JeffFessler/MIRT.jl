# linear-map/getindex.jl
#
# Provides getindex() capabilities like A[:,j] for LinearMap objects.
#
# Currently this provides only a partial set of the possible ways
# one can use indexing for a matrix, because currently
# LinearMap is not a subtype of an AbstractMatrix.
#
# If LinearMap were a subtype of an AbstractMatrix, then all possible
# indexing will be supported by Base.getindex, albeit very likely
# by quite inefficient iterators.
#
# These capabilities are provided as a user convenience.
# The user must recognize that LinearMap objects are not designed
# for efficient element-wise indexing.
#
# 2018-01-19, Jeff Fessler, University of Michigan

using LinearMaps
using Test

# A[end]
function Base.lastindex(A::LinearMap)
    return prod(size(A))
end

# A[?,end] and A[end,?]
function Base.lastindex(A::LinearMap, d::Integer)
    return size(A, d)
end


# A[i,j]
function Base.getindex(A::LinearMap, i::Integer, j::Integer)
	e = zeros(size(A,2)); e[j] = 1
	tmp = A * e
	return tmp[i]
end

# A[k]
function Base.getindex(A::LinearMap, k::Integer)
	c = CartesianIndices(size(A))[k] # is there a more elegant way?
	return A[c[1], c[2]]
#	return A[c[:]] # fails
#	return A[c...] # fails
#	return A[Tuple(c)...] # works
#	return getindex(A, c[1], c[2]) # works
end

# A[:,j]
# it is crucial to provide this function rather than to inherit from
# Base.getindex(A::AbstractArray, ::Colon, ::Integer)
# because Base.getindex does this by iterating (I think).
function Base.getindex(A::LinearMap, ::Colon, j::Integer)
	e = zeros(size(A,2)); e[j] = 1
	return A * e
end

# A[i,:]
function Base.getindex(A::LinearMap, i::Integer, ::Colon)
	# in Julia: A[i,:] = A'[:,i] for real matrix A else need conjugate
	return eltype(A) <: Complex ? conj.(A'[:,i]) : A'[:,i]
end

# A[:,j:k]
# this one is also important for efficiency
function Base.getindex(A::LinearMap, ::Colon, ur::UnitRange)
	return hcat([A[:,j] for j in ur]...)
end

# A[i:k,:]
Base.getindex(A::LinearMap, ur::UnitRange, ::Colon) = A'[:,ur]'

# A[:,:] = Matrix(A)
Base.getindex(A::LinearMap, ::Colon, ::Colon) = Matrix(A)

# A[i:k,j:l]
# this one is inefficient so could be inherited from Base?
Base.getindex(A::LinearMap, r1::UnitRange, r2::UnitRange) = A[:,r2][r1,:]

# A[???]
# informative error message in case we have overlooked any types
#=
function Base.getindex(A::LinearMap, kw...)
	@show kw
	for arg in kw
		@show typeof(arg)
	end
	throw("unsupported indexing type")
end
=#


# tests for getindex for LinearMaps

function ir_LinearMap_test_getindex(A::LinearMap)
	B = Matrix(A)
	@test all(size(A) .>= (4,4)) # required by tests
	@test B[1] == A[1]
	@test B[7] == A[7]
	@test B[:,5] == A[:,5]
	@test B[3,:] == A[3,:]
	@test B[1,3] == A[1,3]
	@test B[:,1:3] == A[:,1:3]
	@test B[1:3,:] == A[1:3,:]
	@test B[1:3,2:4] == A[1:3,2:4]
	@test B == A[:,:]
	@test B'[3] == A'[3]
	@test B'[:,4] == A'[:,4]
	@test B'[2,:] == A'[2,:]
	@test B[end] == A[end] # lastindex
	@test B[3,end-1] == A[3,end-1]

	# The following do not work because currently LinearMap is not a
	# subtype of AbstractMatrix.  If it were such a subtype, then LinearMap
	# would inherit general Base.getindex abilities
#=
	# todo later, if/when LinearMap is a subtype of AbstractMatrix
	@test B[[1, 3, 4]] == A[[1, 3, 4]]
	@test B[:, [1, 3, 4]] == A[:, [1, 3, 4]]
	@test B[[1, 3, 4], :] == A[[1, 3, 4], :]
	@test B[4:7] == A[4:7]
=#

	true
end


# tests for cumsum()
# note: the adjoint of cumsum is reverse(cumsum(reverse(y)))
function ir_lm_test_getindex_cumsum( ; N::Integer=5)
	A = LinearMap(cumsum, y -> reverse(cumsum(reverse(y))), N)
	ir_LinearMap_test_getindex(A)
	@test Matrix(A)' == Matrix(A')
#	@show A[:,4]
	true
end
