# Aodwt.jl

using MIRT: Aodwt
using MIRTjim: jim

using LinearMapsAA: LinearMapAM, LinearMapAO
using Test: @test


"""
    Aodwt_show( ; dims::Dims=(64,32), level::Int=3)
show scales
"""
function Aodwt_show( ; dims::Dims = (64, 32), level::Int=3)
	W, scales, mfun = Aodwt(dims, level=level)
	jim(scales)
end


for op in (true, false)
	W,_,_ = Aodwt((8,16) ; level=2, operator=op)
	@test Matrix(W)' == Matrix(W') # check adjoint
	op && isinteractive() && (@show W.wt)
	op && isinteractive() && (@show propertynames(W.wt))
	@test W.level == 2
	@test W.wt.name == "haar"
	@test W isa (op ? LinearMapAO : LinearMapAM)
end

Aodwt_show()
