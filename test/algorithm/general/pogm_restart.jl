# pogm_restart.jl

using MIRT: pogm_restart

using LinearAlgebra: norm, opnorm, I
using Random: seed!
using Test: @test, @test_throws


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

@test_throws ArgumentError pogm_restart(x0, Fcost, f_grad, a2; mom=:bad)

x, _ = pogm_restart(x0, Fcost, f_grad, a2;
	f_mu=0, mom=:fpgm, niter=100, g_prox=g_prox, fun=fun)
pogm_restart(x0, Fcost, f_grad, a2;
	f_mu=0, mom=:pgm, niter=2, g_prox=g_prox)
pogm_restart(x0, Fcost, f_grad, a2;
	f_mu=reg, mom=:fpgm, niter=2, g_prox=g_prox) # with mu
pogm_restart(x0, Fcost, f_grad, a2;
	f_mu=reg, mom=:pgm, niter=2, g_prox=g_prox) # with mu

# todo: @inferred

x, out = pogm_restart(x0, Fcost, f_grad, a2;
	f_mu=0, mom=:pogm, restart=:gr, restart_cutoff=0.,
	bsig=1, niter=100, g_prox=g_prox, fun=fun)
@test x ≈ xh

x, out = pogm_restart(x0, Fcost, f_grad, a2;
	f_mu=reg, mom=:pogm, restart=:gr, restart_cutoff=0., # with mu
	bsig=1, niter=100, g_prox=g_prox, fun=fun)
@test x ≈ xh
