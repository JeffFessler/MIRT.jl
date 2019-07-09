#=
pogm_restart.jl
2017-03-31, Donghwan Kim and Jeff Fessler, University of Michigan
=#

using LinearAlgebra: norm, opnorm
using Random: seed!
using Test: @test, @test_throws


function gr_restart(Fgrad, ynew_yold, restart_cutoff)
	return sum(Float64, real(-Fgrad .* ynew_yold)) <=
		restart_cutoff * norm(Fgrad) * norm(ynew_yold)
end


"""
`x, out = pogm_restart(x0, Fcost, f_grad, f_L;
	f_mu=0, mom=:pogm, restart=:gr, restart_cutoff=0.,
	bsig=1, niter=10, g_prox=(z,c)->z, fun=...)`

Iterative proximal algorithms (PGM=ISTA, FPGM=FISTA, POGM) with restart.

# in
* `x0` initial guess
* `Fcost` function for computing the cost function value ``F(x)``
  - (needed only if `restart == :fr`)
* `f_grad` function for computing the gradient of ``f(x)``
* `f_L` Lipschitz constant of the gradient of ``f(x)``

# option
* `f_mu` strong convexity parameter of ``f(x)``; default 0.
  - if `f_mu > 0`, ``(\\alpha, \\beta_k, \\gamma_k)`` is chosen by Table 1 in [KF18]
* `g_prox` function `g_prox(z,c)` for the proximal operator for ``g(x)``
  - `g_prox(z,c)` computes ``argmin_x 1/2 \\|z-x\\|^2 + c \\, g(x)``
* `mom`	momentum option
  - `:pogm` POGM (fastest); default!
  - `:fpgm` (FISTA), ``\\gamma_k = 0``
  - `:pgm` PGM (ISTA), ``\\beta_k = \\gamma_k = 0``
* `restart` restart option
  - `:gr` gradient restart; default!
  - `:fr` function restart
  - `:none` no restart
* `restart_cutoff` for `:gr` restart if cos(angle) < this; default 0.
* `bsig` gradient "gamma" decrease option (value within [0 1]); default 1
  - see ``\\bar{\\sigma}`` in [KF18]
* `niter` number of iterations; default 10
* `fun` function`(iter, xk, yk, is_restart)` user-defined function evaluated each `iter` with secondary `xk`, primary `yk`, and boolean `is_restart` indicating whether this iteration was a restart

# out
* `x`	final iterate
  - for ISTA / FPGM (FISTA): primary iterate ``y_N``
  - for POGM: secondary iterate ``x_N``, see [KF18]
* `out [fun(0, x0, x0, false), fun(1, x1, y1, is_restart), ...]` array of length `[niter+1]`

Optimization Problem: Nonsmooth Composite Convex Minimization
* ``argmin_x F(x),  F(x) := f(x) + g(x))``
  - ``f(x)`` smooth convex function
  - ``g(x)`` convex function, possibly nonsmooth and "proximal-friendly" [CP11]

# Optimization Algorithms:

Accelerated First-order Algorithms when ``g(x) = 0`` [KF18]
iterate as below for given coefficients ``(\\alpha, \\beta_k, \\gamma_k)``
* For k = 0,1,...
  - ``y_{k+1} = x_k - \\alpha  f'(x_k)`` : gradient update
  - ``x_{k+1} = y_{k+1} + \\beta_k  (y_{k+1} - y_k) + \\gamma_k  (y_{k+1} - x_k)`` : momentum update

Proximal versions of the above for ``g(x) \\neq 0`` are in the below references,
and use the proximal operater
``prox_g(z) = argmin_x {1/2\\|z-x\\|^2 + g(x)}``.

* Proximal Gradient method (PGM or ISTA) - ``\\beta_k = \\gamma_k = 0``. [BT09]
* Fast Proximal Gradient Method (FPGM or FISTA) - ``\\gamma_k = 0``. [BT09]
* Proximal Optimized Gradient Method (POGM) - [THG15]
* FPGM(FISTA) with Restart - [OC15]
* POGM with Restart - [KF18]

# references

* [CP11] P. L. Combettes, J. C. Pesquet,
 "Proximal splitting methods in signal processing,"
 Fixed-Point Algorithms for Inverse Problems in Science and Engineering,
 Springer, Optimization and Its Applications, 2011.
* [KF18] D. Kim, J.A. Fessler,
 "Adaptive restart of the optimized gradient method for convex optimization," 2018
 Arxiv:1703.04641,
 [http://doi.org/10.1007/s10957-018-1287-4]
* [BT09] A. Beck, M. Teboulle:
 "A fast iterative shrinkage-thresholding algorithm for linear inverse problems,"
 SIAM J. Imaging Sci., 2009.
* [THG15] A.B. Taylor, J.M. Hendrickx, F. Glineur,
 "Exact worst-case performance of first-order algorithms
 for composite convex optimization," Arxiv:1512.07516, 2015,
 SIAM J. Opt. 2017
 [http://doi.org/10.1137/16m108104x]

Copyright 2017-3-31, Donghwan Kim and Jeff Fessler, University of Michigan
2018-08-13 Julia 0.7.0
2019-02-24 interface redesign
"""
function pogm_restart(x0, Fcost::Function, f_grad::Function, f_L::Real;
		f_mu::Real = 0.,
		mom::Symbol = :pogm, # :ogm :gm
		restart::Symbol = :gr, # :fr :none
		restart_cutoff::Real = 0.,
		bsig::Real = 1,
		niter::Integer = 10,
		g_prox::Function = (z, c::Real) -> z,
		fun::Function = (iter::Integer, xk, yk, is_restart::Bool) -> 0)

	!in(mom, (:pgm, :fpgm, :pogm)) && throw(ArgumentError("mom $mom"))
	!in(restart, (:none, :gr, :fr)) && throw(ArgumentError("restart $restart"))
	f_L < 0 && throw(ArgumentError("f_L=$f_L < 0"))
	f_mu < 0 && throw(ArgumentError("f_mu=$f_mu < 0"))
	bsig < 0 && throw(ArgumentError("bsig=$bsig < 0"))
	!((-1 < restart_cutoff) && (restart_cutoff < 1)) &&
		throw(ArgumentError("restart_cutoff=$restart_cutoff"))

	L = f_L
	mu = f_mu
	q = mu/L

	# initialize parameters
	told = 1
	sig = 1
	zetaold = 1 # dummy

	# initialize x
	xold = x0
	yold = x0
	uold = x0
	zold = x0
	Fcostold = Fcost(x0)
	Fgradold = zeros(size(x0)) # dummy

	# save initial
	out = Array{Any}(undef, niter+1)
	out[1] = fun(0, x0, x0, false)

	xnew = []
	ynew = []

# iterations
for iter=1:niter

	# proximal gradient (PGM) update
	if mom == :pgm && mu != 0
		alpha = 2. / (L+mu)
	else
		alpha = 1. / L
	end

	fgrad = f_grad(xold)

	is_restart = false

	if mom == :pgm || mom == :fpgm
		ynew = g_prox(xold - alpha * fgrad, alpha) # standard PG update
		Fgrad = -(1. / alpha) * (ynew - xold) # standard composite gradient mapping
		Fcostnew = Fcost(ynew)

		# restart condition
		if restart != :none
			# function/gradient restart
			if ((restart == :fr && Fcostnew > Fcostold)
			|| (restart == :gr && gr_restart(Fgrad, ynew-yold, restart_cutoff)))
				told = 1
				is_restart = true
			end
			Fcostold = Fcostnew
		end

	elseif mom == :pogm # POGM
		# gradient update for POGM [see KF18]
		unew = xold - alpha * fgrad
		# restart + "gamma" decrease conditions checked later for POGM,
		# unlike PGM, FPGM above

#	else
#		throw("bad mom $mom")
	end

	# momentum coefficient "beta"
	if mom == :fpgm && mu != 0 # known μ > 0
		beta = (1 - sqrt(q)) / (1 + sqrt(q))
	elseif mom == :pogm && mu != 0
		beta = (2 + q - sqrt(q^2+8*q))^2 / 4. / (1-q)
	# for "mu" = 0 or for unknown "mu"
	elseif mom != :pgm
		if mom == :pogm && iter == niter # && restart == 0
			tnew = 0.5 * (1 + sqrt(1 + 8 * told^2))
		else
			tnew = 0.5 * (1 + sqrt(1 + 4 * told^2))
		end

		beta = (told - 1) / tnew
	end

	# momentum update
	if mom == :pgm
		xnew = ynew
	elseif mom == :fpgm
		xnew = ynew + beta * (ynew - yold)
	elseif mom == :pogm # see [KF18]
		# momentum coefficient "gamma"
		if mu != 0
			gamma = (2 + q - sqrt(q^2+8*q)) / 2.
		else
			gamma = sig * told / tnew
		end

		znew = (unew + beta * (unew - uold) + gamma * (unew - xold)
				- beta * alpha / zetaold * (xold - zold))
		zetanew = alpha * (1 + beta + gamma)
		xnew = g_prox(znew, zetanew) # non-standard PG update for POGM

		# non-standard composite gradient mapping for POGM:
		Fgrad = fgrad - 1/zetanew * (xnew - znew)
		ynew = xold - alpha * Fgrad
		Fcostnew = Fcost(xnew)

		# restart + "gamma" decrease conditions for POGM
		if restart != :none
			# function/gradient restart
			if ((restart == :fr && Fcostnew > Fcostold)
			|| (restart == :gr && gr_restart(Fgrad, ynew-yold, restart_cutoff)))
				tnew = 1
				sig = 1
				is_restart = true

			# gradient "gamma" decrease
			elseif sum(Float64, real(Fgrad .* Fgradold)) < 0
				sig = bsig * sig
			end

			Fcostold = Fcostnew
			Fgradold = Fgrad
		end

		uold = unew
		zold = znew
		zetaold = zetanew
	end

	out[iter+1] = fun(iter, xnew, ynew, is_restart) # save

	xold = xnew
	yold = ynew

	if mom != :pgm && mu == 0
		told = tnew
	end
end # for iter

	return ((mom == :pogm) ? xnew : ynew), out
end # pogm_restart()


"""
`pogm_restart_test()`
self test
"""
function pogm_restart_test()
	# Tikhonov regularized LS problem
	seed!(0); M = 30; N = 6; A = randn(M,N); y = randn(M)
	a2 = opnorm(A)^2
	reg = 0.1 * a2
	xh = (A'A + reg*I) \ A'y
	Fcost = (x) -> 1/2 * norm(A * x - y)^2
	cost = (x) -> Fcost(x) + reg/2 * norm(x)^2
	fun = (iter, x, y, is_restart) ->
		(cost(x) - cost(xh), norm(x - xh) / norm(xh), time())
	f_grad = (x) -> A' * (A * x - y)
	g_prox = (z,c) -> z / (1 + reg * c) # proximal operator for 2-norm
	grad = (x) -> A' * (A * x - y) + reg * x
	x0 = A \ y

	x, out = pogm_restart(x0, Fcost, f_grad, a2;
		f_mu=0, mom=:pogm, restart=:gr, restart_cutoff=0.,
		bsig=1, niter=100, g_prox=g_prox, fun=fun)
	@test isapprox(x, xh)
	x, out = pogm_restart(x0, Fcost, f_grad, a2;
		f_mu=reg, mom=:pogm, restart=:gr, restart_cutoff=0., # with mu
		bsig=1, niter=100, g_prox=g_prox, fun=fun)
	@test isapprox(x, xh)
	x, _ = pogm_restart(x0, Fcost, f_grad, a2;
		f_mu=0, mom=:fpgm, niter=100, g_prox=g_prox, fun=fun)
	pogm_restart(x0, Fcost, f_grad, a2;
		f_mu=0, mom=:pgm, niter=2, g_prox=g_prox)
	pogm_restart(x0, Fcost, f_grad, a2;
		f_mu=reg, mom=:fpgm, niter=2, g_prox=g_prox) # with mu
	pogm_restart(x0, Fcost, f_grad, a2;
		f_mu=reg, mom=:pgm, niter=2, g_prox=g_prox) # with mu
	@test_throws ArgumentError pogm_restart(x0, Fcost, f_grad, a2; mom=:bad)
	true
end


"""
`pogm_restart(:test)` self test
"""
function pogm_restart(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	@test pogm_restart_test()
	true
end
