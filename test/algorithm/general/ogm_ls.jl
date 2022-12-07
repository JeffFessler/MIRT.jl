# test/algorithm/general/ogm_ls.jl

using MIRT: ogm_ls

using LinearAlgebra: norm, opnorm, I
using Random: seed!
using Plots: plot, plot!, scatter!, gui
using LaTeXStrings
using Test: @test, @test_throws, @inferred


seed!(0); M = 30; N = 6; A = randn(M,N); y = randn(M)
a2 = opnorm(A)^2
reg = 0.1 * a2
xh = (A'A + reg*I) \ A'y # regularized LS solution
cost = (x) -> 1/2 * norm(A * x - y)^2 + reg/2 * norm(x)^2
fun = (x,iter) -> (cost(x) - cost(xh), norm(x - xh) / norm(xh), time())
grad1 = (x) -> A' * (A * x - y) + reg * x
curv1 = (x) -> a2 + reg

B = [A, I] # matrix blocks
gradf = [u -> u - y, v -> reg * v] # f functions gradients
curvf = [v -> 1, v -> reg]

niter = 90
# (x1, out1) = @inferred ogm_ls(   grad1, curv1, zeros(N); niter, fun) # todo: fails
x1, out1 = ogm_ls(   grad1, curv1, zeros(N); niter, fun)
x2, out2 = ogm_ls(B, gradf, curvf, zeros(N); niter, fun)

#@test x1 ≈ xh # no, converges too slowly
@test x2 ≈ xh

if isinteractive()
    lf = x -> log10(max(x,1e-17))
    costk = out -> lf.([out[k][1] for k=1:niter+1])
    errk = out -> lf.([out[k][2] for k=1:niter+1])
    allk = out -> (costk(out), errk(out))
    cost1, err1 = allk(out1)
    cost2, err2 = allk(out2)

    k = 0:niter
    p1 = plot(xlabel="k", ylabel=L"\log(\Psi(x_k) - \Psi(x_*))")
    scatter!(k, cost1, color=:blue, label="cost1")
    scatter!(k, cost2, color=:red, marker=:x, label="cost2")

    p2 = plot(xlabel="k", ylabel=L"\log(\|x_k - x_*\|/\|x_*\|)")
    scatter!(k, err1, color=:blue, label="NRMSD1")
    scatter!(k, err2, color=:red, marker=:x, label="NRMSD2")
    plot(p1, p2)
    gui()
end
