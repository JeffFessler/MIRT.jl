# ncg.jl

using MIRT: ncg

using LinearAlgebra: I, norm, opnorm
using Plots: Plot
using Random: seed!
using Plots: plot, plot!, scatter!, gui
using LaTeXStrings
using Test: @test, @test_throws, @inferred
using LinearMapsAA


odim = (10,4); idim = (5,2)
M = prod(odim)
N = prod(idim)
seed!(0)
A = randn(M,N)
a2 = opnorm(A)^2
y = randn(M)
reg = 0.1 * a2
xh = (A'A + reg*I) \ A'y
if true # stress test with AO
	A = LinearMapAA(A ; odim=odim, idim=idim)
	xh = reshape(xh, idim)
	y = reshape(y, odim)
end
cost = (x) -> 1/2 * norm(A * x - y)^2 + reg/2 * norm(x)^2
fun = (x,iter) -> (cost(x) - cost(xh), norm(x - xh) / norm(xh), time())
grad1 = (x) -> A' * (A * x - y) + reg * x
curv1 = (x) -> a2 + reg

B = [A, I] # matrix blocks
gradf = [u -> u - y, v -> reg * v] # f functions gradients
curvf = [v -> 1, v -> reg]

#	todo-i: these next two fail:
#	@inferred ncg([ones(3,3)], [x -> x], [v -> 1], zeros(Float32, 3), niter=2)
#	@inferred ncg(x -> x, v -> 1, zeros(Float32, 3), niter=2)
#	error("todo")

niter = 40
x0 = zeros(size(xh))
x1, out1 = ncg(   grad1, curv1, x0, niter=niter, fun=fun)
x2, out2 = ncg(B, gradf, curvf, x0, niter=niter, fun=fun)
@test x1 ≈ xh
@test x2 ≈ xh

@test_throws ArgumentError ncg(grad1, curv1, x0, betahow=:test)

lf = x -> log10(max(x,1e-17))
costk = out -> lf.([out[k][1] for k=1:niter+1])
errk = out -> lf.([out[k][2] for k=1:niter+1])
allk = out -> (costk(out), errk(out))
cost1, err1 = allk(out1)
cost2, err2 = allk(out2)

if isinteractive()
	k = 0:niter
	plot(xlabel="k", ylabel=L"\log(\Psi(x_k) - \Psi(x_*))")
	scatter!(k, cost1, color=:blue, label="cost1")
	scatter!(k, cost2, color=:red, marker=:x, label="cost2")
	p1 = plot!()

	plot(xlabel="k", ylabel=L"\log(\|x_k - x_*\|/\|x_*\|)")
	scatter!(k, err1, color=:blue, label="NRMSD1")
	scatter!(k, err2, color=:red, marker=:x, label="NRMSD2")
	p2 = plot!()
	plot(p1, p2)
	gui()
end


#=
todo
f = x -> x + one(eltype(x))
f = x -> eltype(x).(x)

@inferred f(3)
@inferred f(3.0f0)

function g(x0::AbstractVector{<:Number} ; f::Function = x -> 2*x)
	x = copy(x0)
#	return x .+ ones(eltype(x0), size(x0))
#	return 2*x
	return f(x)
end

g(1:3)
@inferred g(1:3)
@code_warntype g(1:3)

#	the following test returns nothing!?
#	@code_warntype ncg([x -> 2x], [v -> 1], zeros(Float32, 3), niter=2)
=#
