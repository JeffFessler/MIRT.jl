# ogm_ls.jl
# OGM with a MM line search
# 2019-03-16, Jeff Fessler, University of Michigan

using LinearAlgebra: I, norm

"""
`(x,out) = ogm_ls(B, gradf, curvf, x0; niter=?, ninner=?, fun=?)`

OGM with a line search; Drori&Taylor @arxiv 1803.05676;
to minimize a general "inverse problem" cost function of the form
``\\Psi(x) = \\sum_{j=1}^J f_j(B_j x)``
where each function ``f_j(v)`` has a quadratic majorizer of the form
`` q_j(v;u) = f_j(u) + \\nabla f_j(u) (v - u) + 1/2 \\|v - \\u|^2_{C_j(u)} ``
where ``C_j(u)`` is diagonal matrix of curvatures.
(It suffices for each ``f_j`` to have a Lipschitz smooth gradient.)

This OGM method uses a majorize-minimize (MM) line search.

in
* `B`		array of ``J`` blocks ``B_1,...,B_J``
* `gradf`	array of ``J`` functions return gradients of ``f_1,...,f_J``
* `curvf`	array of ``J`` functions `z -> curv(z)` that return a scalar
		or a vector of curvature values for each element of ``z``
* `x0`	initial guess; need `length(x) == size(B[j],2)` for ``j=1...J``

option
* `niter`	# number of outer iterations; default 50
* `ninner`	# number of inner iterations of MM line search; default 5
* `fun`		User-defined function to be evaluated with two arguments (x,iter).
			It is evaluated at (x0,0) and then after each iteration.

output
* `x`		final iterate
* `out`		`[niter+1] (fun(x0,0), fun(x1,1), ..., fun(x_niter,niter))`
	(all 0 by default). This is an array of length `niter+1`
"""
function ogm_ls(
	B::AbstractVector{<:Any},
	gradf::AbstractVector{<:Function},
	curvf::AbstractVector{<:Function},
	x0::AbstractVector{<:Number};
	niter::Integer=50,
	ninner::Integer=5,
	fun::Function = (x,iter) -> 0)

out = Array{Any}(undef, niter+1)
out[1] = fun(x0, 0)

J = length(B)

x = x0
dir = []
grad_old = []
grad_new = []
grad_sum = zeros(size(x0))
ti = 1
thetai = 1

B0 = [B[j] * x for j=1:J]
Bx = copy(B0)
By = copy(B0)
grad = (Bx) -> sum([B[j]' * gradf[j](Bx[j]) for j=1:J])

for iter = 1:niter
	grad_new = grad(Bx) # gradient of x_{iter-1}
	grad_sum += ti * grad_new # sum_{j=0}^{iter-1} t_j * gradient_j

	thetai = (1 + sqrt(8*ti^2 + 1)) / 2 # theta_{i+1}
	ti = (1 + sqrt(4*ti^2 + 1)) / 2 # t_{i+1}

	tt = (iter < niter) ? ti : thetai # use theta_i factor for last iteration
	yi = (1 - 1/tt) * x + (1/tt) * x0

	for j=1:J # update Bj * yi
		By[j] = (1 - 1/tt) * Bx[j] + (1/tt) * B0[j]
	end

	dir = -(1 - 1/tt) * grad_new - (2/tt) * grad_sum # -d_i

	# MM-based line search for step size alpha
	# using h(a) = sum_j f_j(By_j + a * Bd_j)
	Bd = [B[j] * dir for j=1:J]

	alf = 0
	for ii=1:ninner
		derh = 0 # derivative of h(a)
		curv = 0
		for j=1:J
			tmp = By[j] + alf * Bd[j]
			derh += real(Bd[j]' * gradf[j](tmp))
			curv += sum(curvf[j](tmp) .* abs.(Bd[j]).^2)
		end
		curv < 0 && throw("curv < 0")
		if curv > 0
			alf = alf - derh / curv
		end
		if alf == 0
			break
		end
	end

#	# derivative of h(a) = cost(x + a * dir) where \alpha is real
#	dh = alf -> real(sum([Bd[j]' * gradf[j](By[j] + alf * Bd[j]) for j=1:J]))
#	Ldh = sum([Lgf[j] * norm(Bd[j])^2 for j=1:J]) # Lipschitz constant for dh
#	(alf, ) = gd(dh, Ldh, 0, niter=ninner) # GD-based line search
# todo

	x = yi + alf * dir

	if iter < niter
		for j=1:J # update Bj * x
			Bx[j] = By[j] + alf * Bd[j]
		end
	end

#	for j=1:J # recursive update Bj * yi ???
#		By[j] = (1 - 1/ti) * (By[j] + alf * Bd[j]) + (1/ti) * B0[j]
#	end

	out[iter+1] = fun(x, iter)
end

return x, out
end


"""
`(x,out) = ogm_ls(grad, curv, x0, ...)`

special case of `ogm_ls` (OGM with line search) for minimizing a cost function
whose gradient is `grad(x)`
and that has a quadratic majorizer with diagonal Hessian given by
`curv(x)`.
Typically `curv = (x) -> L` where `L` is the Lipschitz constant of `grad`
"""
function ogm_ls(
	grad::Function,
	curv::Function,
	x0::AbstractVector{<:Number};
	kwargs...)

	return ogm_ls([I], [grad], [curv], x0; kwargs...)
end


using LinearAlgebra: norm, opnorm, I
using Random: seed!
using Plots
using LaTeXStrings

function ogm_ls_test()
	seed!(0); M = 30; N = 6; A = randn(M,N); y = randn(M)
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

	niter = 90
	x1, out1 = ogm_ls(   grad1, curv1, zeros(N), niter=niter, fun=fun)
	x2, out2 = ogm_ls(B, gradf, curvf, zeros(N), niter=niter, fun=fun)

	lf = x -> log10(max(x,1e-17))
	costk = out -> lf.([out[k][1] for k=1:niter+1])
	errk = out -> lf.([out[k][2] for k=1:niter+1])
	allk = out -> (costk(out), errk(out))
	cost1, err1 = allk(out1)
	cost2, err2 = allk(out2)

#	!isapprox(x1, xh) && throw("bug: x1 vs xh") # no, converges too slow
	!isapprox(x2, xh) && throw("bug: x2 vs xh")

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
`ogm_ls(:test)`

run test
"""
function ogm_ls(test::Symbol)
	test != :test && throw("test")
	ogm_ls_test()
	true
end
