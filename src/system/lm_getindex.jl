# lm_getindex.jl
#
# Provides getindex() capabilities like A[:,j] for LinearMap objects
# Currently this provides only a partial set of the possible ways
# one can use indexing for a matrix, because currently
# LinearMap is not a subtype of an AbstractMatrix.
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
	return isreal(A) ? A'[:,i] : conj.(A'[:,i])
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
if false
function Base.getindex(A::LinearMap, kw...)
	@show kw
	for arg in kw
		@show typeof(arg)
	end
	error("unsupported indexing type")
end
end


# tests for getindex for LinearMaps

function ir_LinearMap_test_getindex(A::LinearMap)
	B = Matrix(A)
	@assert all(size(A) .>= (4,4)) # required by tests
	@assert B[1] == A[1]
	@assert B[7] == A[7]
	@assert B[:,5] == A[:,5]
	@assert B[3,:] == A[3,:]
	@assert B[1,3] == A[1,3]
	@assert B[:,1:3] == A[:,1:3]
	@assert B[1:3,:] == A[1:3,:]
	@assert B[1:3,2:4] == A[1:3,2:4]
	@assert B == A[:,:]
	@assert B'[3] == A'[3]
	@assert B'[:,4] == A'[:,4]
	@assert B'[2,:] == A'[2,:]

	# The following do not work because currently LinearMap is not a
	# subtype of AbstractMatrix.  If it were such a subtype, then LinearMap
	# would inherit general Base.getindex abilities
	if false # todo later, if/when LinearMap is a subtype of AbstractMatrix
		@assert B[[1, 3, 4]] == A[[1, 3, 4]]
		@assert B[:, [1, 3, 4]] == A[:, [1, 3, 4]]
		@assert B[[1, 3, 4], :] == A[[1, 3, 4], :]
		@assert B[4:7] == A[4:7]
	end
	true
end


# tests for cumsum()
# note: the adjoint of cumsum is reverse(cumsum(reverse(y)))
function ir_lm_test_getindex_cumsum(;N::Integer=5)
	A = LinearMap(cumsum, y -> reverse(cumsum(reverse(y))), N)
	ir_LinearMap_test_getindex(A)
	@assert Matrix(A)' == Matrix(A')
#	@show A[:,4]
	true
end
