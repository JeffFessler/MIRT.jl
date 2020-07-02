#=
Aodwt.jl
2019-02-23 Jeff Fessler, University of Michigan
=#

export Aodwt

using Plots
using LinearMapsAA: LinearMapAA, LinearMapAM, LinearMapAO
using Wavelets: dwt!, idwt!, wavelet, WT
using Test: @test
#using MIRT: jim


"""
    A, levels, mfun = Aodwt(dims ; level::Int=3, wt=wavelet(WT.haar))

Create orthogonal discrete wavelet transform (ODWT) `LinearMapAA`

in
- `dims::Dims` tuple of dimensions

option
- `level::Int` # of levels; default 3
- `wt` wavelet transform type (see `Wavelets` package); default Haar
- `operator::Bool=true` default to `LinearMapAO`

out
- `A` a `LinearMapAX` object
- `scales` array of size `dims` showing the scale of each coefficient
which is useful when imposing scale-dependent regularization
- `mfun` convenience function for A*X when X is a Matrix or Array (not vector)

2019-02-23 Jeff Fessler, University of Michigan
"""
function Aodwt(
	dims::Dims ;
	T::DataType = Float32,
	level::Int = 3,
	wt = wavelet(WT.haar),
	operator::Bool = true, # !
)

	function mfunA(lev)
		# todo: avoiding convert involves `dwt!` limitation; see:
		# https://github.com/JuliaDSP/Wavelets.jl/issues/53
		# https://github.com/JuliaDSP/Wavelets.jl/pull/54
	#	forw!(y,x) = dwt!(y, convert(AbstractArray{T},x), wt, level)
	#	back!(x,y) = idwt!(x, convert(AbstractArray{T},y), wt, level)
		forw!(y,x) = dwt!(y, x, wt, lev)
		back!(x,y) = idwt!(x, y, wt, lev)

		if operator
			mfun = (A, x) -> A * x
			return mfun, LinearMapAA(forw!, back!,
				(prod(dims), prod(dims)) ;
				prop = (wt=wt, level=lev), T=T,
				operator = true, idim=dims, odim=dims,
			)
		else
			mfun = (A, x) -> reshape(A * vec(x), dims)
			return mfun, LinearMapAA(
				(y,x) -> vec(forw!(reshape(y, dims), reshape(x, dims))),
				(x,y) -> vec(back!(reshape(x, dims), reshape(y, dims))),
			#	x -> vec(dwt(reshape(x, dims), wt, level)),
			#	y -> vec(idwt(reshape(y, dims), wt, level)),
				(prod(dims), prod(dims)) ;
				prop = (wt=wt, level=lev), T=T,
			)
		end
	end

	mfun, A = mfunA(level)

	scales = zeros(dims)
	for il=1:level
		_,Al = mfunA(il)
		tmp = mfun(Al, ones(dims)) .== 0
		scales += il * (tmp .& (scales .== 0))
	end

	return A, scales, mfun
end


function Aodwt()
	@doc Aodwt
end


"""
    Aodwt_show( ; dims::Dims=(64,32), level::Int=3)
show scales
"""
function Aodwt_show( ; dims::Dims = (64, 32), level::Int=3)
	W, scales, mfun = Aodwt(dims, level=level)
	jim(scales)
end


"""
    Aodwt(:test)
self test

`Aodwt(:show)`
visualize
"""
function Aodwt(test::Symbol)
	if test === :show
		return Aodwt_show()
	end

	test != :test && throw(ArgumentError("test $test"))
	Aodwt()

	for op in (true, false)
		W,_,_ = Aodwt((8,16) ; level=2, operator=op)
		@test Matrix(W)' == Matrix(W') # check adjoint
		op && isinteractive() && (@show W.wt)
		op && isinteractive() && (@show propertynames(W.wt))
		@test W.level == 2
		@test W.wt.name == "haar"
		@test W isa (op ? LinearMapAO : LinearMapAM)
	end

	Aodwt(:show)
	true
end
