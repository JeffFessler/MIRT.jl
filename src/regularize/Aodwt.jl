#=
Aodwt.jl
2019-02-23 Jeff Fessler, University of Michigan
=#

export Aodwt

using Plots
using LinearMapsAA: LinearMapAA
using Wavelets: dwt, idwt, wavelet, WT
#using MIRT: jim


"""
`A, levels, mfun = Aodwt(dims ; level::Int=3, wt=wavelet(WT.haar))`

create orthogonal discrete wavelet transform (ODWT) `LinearMapAA`

in
- `dims::Dims` tuple of dimensions

option
- `level::Int` # of levels; default 3
- `wt` wavelet transform type (see `Wavelets` package); default Haar

out
- `A` a `LinearMapAA` object
- `scales` array of size `dims` showing the scale of each coefficient
which is useful when imposing scale-dependent regularization
- `mfun` convenience function for A*X when X is a Matrix or Array (not vector)

2019-02-23 Jeff Fessler, University of Michigan
"""
function Aodwt(dims::Dims ; level::Int=3, wt=wavelet(WT.haar))

	Afun = (level) -> LinearMapAA(
		x -> dwt(reshape(x, dims), wt, level)[:],
		y -> idwt(reshape(y, dims), wt, level)[:],
		(prod(dims), prod(dims)), (wt=wt, level=level))

	A = Afun(level)

	mfun = (A, x) -> reshape(A * x[:], dims)

	scales = zeros(dims)
	for il=1:level
		Al = Afun(il)
		tmp = mfun(Al, ones(dims)) .== 0
		scales += il * (tmp .& (scales .== 0))
	end

	return A, scales, mfun
end


function Aodwt()
	@doc Aodwt
end


"""
`Aodwt_show( ; M::Int=32, N::Int=64)`
show scales
"""
function Aodwt_show(; dims::Dims = (32, 64))
	W, scales, mfun = Aodwt(dims)
	jim(scales)
end


"""
`Aodwt(:test)`
self test

`Aodwt(:show)`
visualize
"""
function Aodwt(test::Symbol)
	if test == :show
		return Aodwt_show()
	end

	test != :test && throw(ArgumentError("test $test"))
	Aodwt()

	W,_,_ = Aodwt((8,16), level=2)
	@test Matrix(W)' == Matrix(W') # check adjoint
#	@show W.wt
#	@show propertynames(W.wt)
	@test W.level == 2
	@test W.wt.name == "haar"

	Aodwt(:show)
	true
end
