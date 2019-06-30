# ncg.jl
# Nonlinear CG optimization
# 2019-03-16, Jeff Fessler, University of Michigan

using LinearAlgebra: I, norm
using Test: @test, @test_throws
using Plots: Plot


"""
`(x,out) = ncg(B, gradf, curvf, x0; ...)`

Nonlinear preconditioned conjugate gradient algorithm
to minimize a general "inverse problem" cost function of the form
``\\Psi(x) = \\sum_{j=1}^J f_j(B_j x)``
where each function ``f_j(v)`` has a quadratic majorizer of the form
`` q_j(v;u) = f_j(u) + \\nabla f_j(u) (v - u) + 1/2 \\|v - \\u|^2_{C_j(u)} ``
where ``C_j(u)`` is diagonal matrix of curvatures.
(It suffices for each ``f_j`` to have a Lipschitz smooth gradient.)

This CG method uses a majorize-minimize (MM) line search.

in
* `B`		array of ``J`` blocks ``B_1,...,B_J``
* `gradf`	array of ``J`` functions return gradients of ``f_1,...,f_J``
* `curvf`	array of ``J`` functions `z -> curv(z)` that return a scalar
		or a vector of curvature values for each element of ``z``
* `x0`	initial guess; need `length(x) == size(B[j],2)` for ``j=1...J``

option
* `niter`	# number of outer iterations; default 50
* `ninner`	# number of inner iterations of MM line search; default 5
* `P`		# preconditioner; default `I`
* `betahow`	"beta" method for the search direction; default `:dai_yuan`
* `fun`		User-defined function to be evaluated with two arguments (x,iter).
			It is evaluated at (x0,0) and then after each iteration.

output
* `x`		final iterate
* `out`		`[niter+1] (fun(x0,0), fun(x1,1), ..., fun(x_niter,niter))`
	(all 0 by default). This is an array of length `niter+1`
"""
function ncg(
	B::AbstractVector{<:Any},
	gradf::AbstractVector{<:Function},
	curvf::AbstractVector{<:Function},
	x0::AbstractVector{<:Number};
	niter::Integer=50,
	ninner::Integer=5,
	P=I,
	betahow::Symbol=:dai_yuan,
	fun::Function = (x,iter) -> 0)

out = Array{Any}(undef, niter+1)
out[1] = fun(x0, 0)

J = length(B)

x = x0
dir = []
grad_old = []
grad_new = []

Bx = [B[j] * x for j=1:J] # u_j in course notes
grad = (Bx) -> sum([B[j]' * gradf[j](Bx[j]) for j=1:J])

for iter = 1:niter
	grad_new = grad(Bx) # gradient
	npgrad = -(P * grad_new)
	if iter == 1
		dir = npgrad
	else
		if betahow == :dai_yuan
			denom =	(grad_new - grad_old)' * dir
			if denom == 0
				betaval = 0
			else
				betaval = grad_new' * (P * grad_new) / denom
			end
		else
			throw(ArgumentError("unknown beta choice: $betahow"))
		end
		dir = npgrad + betaval * dir # search direction
	end
	grad_old = grad_new

	# MM-based line search for step size alpha
	# using h(a) = sum_j f_j(uj + a vj)
	Bd = [B[j] * dir for j=1:J] # v_j in course notes

	alf = 0
	for ii=1:ninner
	#	derh = alf -> sum([Bd[j]' * gradf[j](Bx[j] + alf * Bd[j]) for j=1:J])
		derh = 0 # derivative of h(a)
		curv = 0
		for j=1:J
			tmp = Bx[j] + alf * Bd[j]
			derh += real(Bd[j]' * gradf[j](tmp))
			curv += sum(curvf[j](tmp) .* abs.(Bd[j]).^2)
		end
		curv < 0 && throw("bug: curv < 0")
		if curv > 0
			alf = alf - derh / curv
		end
		if alf == 0
			break
		end
	end

	x += alf * dir
	for j=1:J # update Bj * x
		Bx[j] += alf * Bd[j]
	end
	out[iter+1] = fun(x, iter)
end

return x, out
end


"""
`(x,out) = ncg(grad, curv, x0, ...)`

special case of `ncg` (nonlinear CG) for minimizing a cost function
whose gradient is `grad(x)`
and that has a quadratic majorizer with diagonal Hessian given by
`curv(x)`.
Typically `curv = (x) -> L` where `L` is the Lipschitz constant of `grad`
"""
function ncg(
	grad::Function,
	curv::Function,
	x0::AbstractVector{<:Number};
	kwargs...)

	return ncg([I], [grad], [curv], x0; kwargs...)
end


using LinearAlgebra: norm, opnorm, I
using Random: seed!
using Plots
using LaTeXStrings

function ncg_test()
	seed!(0); M = 40; N = 10; A = randn(M,N); y = randn(M)
	a2 = opnorm(A)^2
	reg = 0.1 * a2
	xh = (A'A + reg*I) \ A'y
	cost = (x) -> 1/2 * norm(A * x - y)^2 + reg/2 * norm(x)^2
	fun = (x,iter) -> (cost(x) - cost(xh), norm(x - xh) / norm(xh), time())
	grad1 = (x) -> A' * (A * x - y) + reg * x
	curv1 = (x) -> a2 + reg

	B = [A, I] # matrix blocks
	gradf = [u -> u - y, v -> reg * v] # f functions gradients
	curvf = [v -> 1, v -> reg]

	niter = 40
	x1, out1 = ncg(   grad1, curv1, zeros(N), niter=niter, fun=fun)
	x2, out2 = ncg(B, gradf, curvf, zeros(N), niter=niter, fun=fun)
	!isapprox(x1, xh) && throw("bug: x1 vs xh")
	!isapprox(x2, xh) && throw("bug: x2 vs xh")

	@test_throws ArgumentError ncg(grad1, curv1, zeros(N), betahow=:test)

	lf = x -> log10(max(x,1e-17))
	costk = out -> lf.([out[k][1] for k=1:niter+1])
	errk = out -> lf.([out[k][2] for k=1:niter+1])
	allk = out -> (costk(out), errk(out))
	cost1, err1 = allk(out1)
	cost2, err2 = allk(out2)

	k = 0:niter
	plot(xlabel="k", ylabel=L"\log(\Psi(x_k) - \Psi(\hat{x}))")
	scatter!(k, cost1, color=:blue, label="cost1")
	scatter!(k, cost2, color=:red, marker=:x, label="cost2")
	p1 = plot!()

	plot(xlabel="k", ylabel=L"\log(\|x_k - \hat{x}\|/\|\hat{x}\|)")
	scatter!(k, err1, color=:blue, label="NRMSD1")
	scatter!(k, err2, color=:red, marker=:x, label="NRMSD2")
	p2 = plot!()
	plot(p1, p2)
end


"""
`ncg(:test)`

run test
"""
function ncg(test::Symbol)
	test != :test && throw("test")
	@test ncg_test() isa Plots.Plot
	true
end
