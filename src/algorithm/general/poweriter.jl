# poweriter.jl
# 2019-06-06, Jeff Fessler, University of Michigan

using Random: seed!
using LinearAlgebra: opnorm, norm
using Test: @test

"""
`v1,σ1 = poweriter(A; niter=?, ...)`

determine first right singular vector `v1`
and first singular value `σ1` of `A`
by applying power iteration to `A'A`

in
* `A` M × N matrix

option
* `niter` default 200
* `x0` initial guess of `v1`
* `tol` stopping tolerance for s1, default 1e-6
* `chat::Bool` verbose? default false

out
* `v1` `[N]` principal right singular vector
* `σ1` spectral norm of `A`
"""
function poweriter(A; niter=200,
		tol::Real = 1e-6,
		x0::AbstractArray{<:Number} = ones(eltype(A), size(A,2)),
		chat::Bool = true,
	)
	x = copy(x0)
	ratio_old = Inf
	for iter=1:niter
		Ax = A * x
		ratio = norm(Ax) / norm(x)
		if abs(ratio - ratio_old) / ratio < tol
			chat && @info "done at iter $iter"
			break
		end
		ratio_old = ratio
		x = A' * Ax
		x /= norm(x)
	end
	return x, norm(A * x) / norm(x)
end


"""
poweriter(:test)
self test
"""
function poweriter(test::Symbol)
	test != :test && throw("bad symbol $test")
	seed!(0)
	M = 30
	N = 20

	A = randn(M,N) # real
	s0 = opnorm(A)
	chat = false
	_,s1 = poweriter(A; tol=1e-9, chat=chat)
	@test isequal(round(s0,digits=7), round(s1,digits=7))

	A = randn(ComplexF32, M, N) # complex
	s0 = opnorm(A)
	_,s1 = poweriter(A; tol=1e-9, chat=chat)
	@test isequal(round(s0,digits=5), round(s1,digits=5))

	true
end
