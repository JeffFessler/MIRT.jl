#=
hcat.jl
provide hcat support for LinearMaps
Jeff Fessler, University of Michigan
with advice from Daniel Karrasch per
https://github.com/Jutho/LinearMaps.jl/issues/45
=#

export hcat_lm

using LinearMaps
using LinearAlgebra: I, UniformScaling
using Base: hcat
using Test

function hcat_lm(As...) # helper routine for hcat of LinearMaps
	nrow = 0
	for i=1:length(As) # find first non-I to specify number of rows
		if As[i] != I
			nrow = size(As[i],1)
			break
		end
	end
	nrow == 0 && throw(DimensionMismatch("need one non-I")) # [I I] illegal

	for i=1:length(As) # check row consistency
		if As[i] != I
			ndims(As[i]) != 2 && throw(DimensionMismatch("need ndims=2"))
			size(As[i],1) != nrow && 
				throw(DimensionMismatch("need same # of rows, " *
					"got (nrow = $nrow, $i => $(size(As[i],1)))"))
		end
	end
	ncol = map(A -> A == I ? nrow : size(A, 2), As)
	ncol = [0; ncol...] # trick for ranger()

	ranger = i -> (1+sum(ncol[1:i])) : sum(ncol[1:(i+1)]) # ith input range
	f = x -> sum([As[i] * (@view x[ranger(i)]) for i=1:length(As)])
	fc = y -> vcat(map(A -> A'*y, As)...)

	lacks_fc = A -> isa(A, LinearMap) && A.fc == nothing
	if any(map(A -> lacks_fc(A), As)) # any LinearMap missing adjoint?
		return LinearMap(f, nrow, sum(ncol)) # no fc
	end
	fc = y -> vcat(map(A -> A'*y, As)...)
	return LinearMap(f, fc, nrow, sum(ncol))
end


# overload hcat for LinearMap
Base.hcat(A::LinearMap, B::LinearMap) = hcat_lm(A, B)

# [A I] [A I I] [A I B C] etc.
Base.hcat(A::LinearMap, Bs...) = hcat_lm(A, Bs...)

Base.hcat(I::UniformScaling{<:Number}, A::LinearMap, Bs...) =
	hcat_lm(I, A, Bs...)

# this will support a single leading matrix, but not more than one
Base.hcat(A::AbstractMatrix{<:Number}, B::LinearMap, Cs...) =
	hcat_lm(A, B, Cs...)


function hcat_lm_test()
	a = ones(3,2)
	b = zeros(3,4)

	A = LinearMap(x -> a*x, y -> a'*y, 3, 2)
	B = LinearMap(x -> b*x, y -> b'*y, 3, 4)

	C = [A B I] # now returns a new LinearMap :)
	@test size(C) == (3,2+4+3)
	@test C * ones(size(C,2)) == 3. * ones(3)
	@test C' * ones(size(C,1)) == [3.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
	@test Matrix(C)' == Matrix(C')

	D = [I A] # LinearMap
	@test Matrix(D)' == Matrix(D')

	E = [a B] # LinearMap
	@test Matrix(E)' == Matrix(E')

	F = LinearMap(x -> b*x, size(b)...)
	G = [A F] # LinearMap
	@test Matrix(G) == [a b]
	true
end


"""
`hcat_lm(:test)`
self test
"""
function hcat_lm(test::Symbol)
	test != :test && throw("hcat_lm test")
	@test hcat_lm_test()
	true
end
