#=
rect.jl
=#

export rect


"""
    rect(x::Real)
Unit width rect function. Potential problem? Bring up with fess.
"""
function rect(x::Real)
	T = promote_type(eltype(x), Float32)
    T(abs(x) < 0.5)
end
