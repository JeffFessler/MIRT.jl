#=
# [NCG for LASSO](@id 4-ncg-lasso)

Example illustrating the Nonlinear Conjugate Gradient (NCG) method
that in turn uses a
line-search method
based on majorize-minimize (MM) principles
in the Julia package
[`MIRT`](https://github.com/JeffFessler/MIRT.jl).

This page was generated from a single Julia file:
[4-ncg-lasso.jl](@__REPO_ROOT_URL__/4-ncg-lasso.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`4-ncg-lasso.ipynb`](@__NBVIEWER_ROOT_URL__/4-ncg-lasso.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`4-ncg-lasso.ipynb`](@__BINDER_ROOT_URL__/4-ncg-lasso.ipynb).


# ### Setup

# Packages needed here.

using Plots; default(markerstrokecolor = :auto, label="")
using MIRTjim: prompt
#using MIRT: line_search_mm, LineSearchMMWork # todo
using MIRT: NCG, ncg
# using LineSearches: BackTracking, HagerZhang, MoreThuente
# using Optim: cg? todo
using LinearAlgebra: norm, dot, I, mul!
using Random: seed!; seed!(0)
using SparseArrays: sprandn
using BenchmarkTools: @btime, @benchmark
using InteractiveUtils: versioninfo

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() && prompt(:prompt);


#=
# Theory

Many methods for solving inverse problems
involve optimization problems
of the form
```math
\hat{x} = \arg\min_{x ∈ \mathbb{F}^N} f(x)
,\qquad
f(x) = \sum_{j=1}^J f_j(B_j x)
```
where ``\mathbb{F}`` denotes the field of real or complex numbers,
matrix ``B_j`` has size ``M_j × N``,
and ``f_j : \mathbb{F}^{M_j} ↦ \mathbb{R}``.

One could apply general-purpose optimization methods here,
like those in
[Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl),
but often we can obtain faster results
by exploiting the specific
(yet still fairly general) structure,
particularly when the problem dimension ``N`` is large.

Here we focus further
on the case where each function ``f_j(⋅)``
has a quadratic majorizer
of the form
```math
f_j(x) ≤ q_j(x,z) = f_j(z) + \text{real}(⟨ ∇f_j(z), x - z ⟩)
+ \frac{1}{2} (x - z)' D_j(z) (x - z),
```
where ``D_j(z)`` is a positive semidefinite matrix
that typically is diagonal.
Often it is a constant times the identity matrix,
e.g.,
a Lipschitz constant for ``∇f_j``,
but often there are
sharper majorizers.
Those majorizers help with efficient line searches.
See 3-ls-mm. # todo

# Smooth LASSO problem

The usual LASSO optimization problem uses the cost function
```math
f(x) = \frac{1}{2} \| A x - y \|_2^2 + β R(x)
,\qquad
R(x) = \| x \|_1 = \sum_{n=1}^N |x_n| = 1' \text{abs.}(x).
```

The 1-norm is just a relaxation of the 0-norm
so here we further "relax" it
by considering the "corner rounded" version
using the Fair potential function
```math
R(x) = \sum_{n=1}^N ψ(x_n) = 1' ψ.(x),
\qquad
ψ(z) = δ^2 |z/δ| - \log(1 + |z/δ|)
```
for a small value of ``δ``.

The derivative of this potential function
is
``\dot{ψ}(z) = z / (1 + |z / δ|)``
and Huber's curvature
``ω_{ψ}(z) = 1 / (1 + |z / δ|)``
provides a suitable majorizer.

This smooth LASSO cost function has the general form above
with ``J=2``,
``B_1 = A``,
``B_2 = I``,
``f_1(u) = \frac{1}{2} \| u - y \|_2^2,``
``f_2(u) = β 1' ψ.(u),``
for which
``∇f_1(u) = u - y,``
``∇f_2(u) = β ψ.(u),``
and
``∇^2 f_1(u) = I,``
``∇^2 f_2(u) \succeq β \, \text{diag}(ω_{ψ}(u)).``
=#


# Fair potential, its derivative and Huber weighting function:
function fair_pot()
    fpot(z,δ) = δ^2 * (abs(z/δ) - log(1 + abs(z/δ)))
    dpot(z,δ) = z / (1 + abs(z/δ))
    wpot(z,δ) = 1 / (1 + abs(z/δ))
    return fpot, dpot, wpot
end;

# Data, cost function and gradients for smooth LASSO problem:
M, N = 1000, 2000
A = sprandn(M, N, 0.3) # sparse system matrix
tmp = sum(abs2, A; dims=1)
tmp = sqrt(sum(tmp) / length(tmp))
A ./= tmp # normalize so that diag(A'A) ≈ 1
#A = I(N)
xtrue = randn(N) .* (rand(N) .< 0.3) # sparse vector
xtrue = sort(xtrue) # since A is random anyway
ytrue = A * xtrue
y = ytrue + 0.0001 * randn(size(ytrue))
snr = 10 * log10(norm(ytrue) / norm(y - ytrue))
# @show snr
β = 2^2
δ = 0.04
fpot, dpot, wpot = Base.Fix2.(fair_pot(), δ)

f(x) = 0.5 * norm(A * x - y)^2 + β * sum(fpot, x)
∇f(x) = A' * (A * x - y) + β * dpot.(x)
x0 = randn(N) # random initial point
x0 = 0 * x0 # zero initial point


# Apply NCG with MM-based line search: simple version.
# The key inputs are the gradient and curvature functions:
B = [A, I]
gradf = [
    u -> u - y, # ∇f₁ for data-fit term
    u -> β * dpot.(u), # ∇f₂ for regularizer
]
curvf = [
    u -> 1, # curvature for data-fit term
    u -> β * wpot.(u), # Huber curvature for regularizer
]
gradf = [
#   let y=y; z -> Iterators.map(-, z, y); end, # z - y
    let y=y, work=similar(y); z -> (@. work = z - y); end, # z - y
    let β=β, dpot=dpot, work=similar(xtrue); z -> (@. work = β * dpot(z)); end, # β * dψ.(z)
]
curvf = [
    u -> 1, # todo only needed for old ncg
#   let β=β, wpot=wpot, work=similar(xtrue); z -> (@. work = β * wpot(z)); end, # β * ωψ.(z)
    let β=β, wpot=wpot; z -> Iterators.map(z -> β * wpot(z), z); end, # β * ωψ.(z)
]

#=
import MIRT
@time MIRT._update!(state);
0.002177 seconds (76 allocations : 253.438 KiB)
0.002279 seconds (56 allocations: 95.500 KiB) # after gradf let
0.002242 seconds (35 allocations: 592 bytes) # after curvf let
=#


niter = 20
ninner = 3
fun = (x, iter) -> f(x) # log this
xcg1, out1 = ncg(B, gradf, curvf, x0; niter, ninner, fun)

pc = plot(0:niter, out1, marker=:circle, label="ncg",
   xlabel="Iteration", ylabel="Cost f(x)")

px = plot(xtrue, label="xtrue")
scatter!(px, xcg1, label="xcg1")
plot!(px, xtrue, label="")
pv = plot(xtrue, xcg1, label="xcg1", aspect_ratio=:equal)
pxv = plot(px, pv)
#gui()
#throw()

#=
todo
out = Vector{Any}(undef, niter+1)
=#

#=
Thanks to NCG and the efficient line search using Huber's curvatures,
the ``f(x_t)`` sequence converges very quickly.

Now compare to `NCG`.
=#

state = NCG(B, gradf, curvf, x0; niter, ninner)

out2 = similar(out1)
out2[1] = f(x0)
for item in state
    out2[state.iter+1] = f(state.x)
end
xcg2 = state.x

# Compare cost functions
# They take slightly different paths, for unknown reasons.
# But the final estimates are quite similar.

plot!(pc, 0:niter, out2, label="NCG", marker=:square)
plot!(pv, xtrue, xcg2, label="xNCG")
scatter!(px, xcg2, label="xNCG")
# scatter(xcg1, xcg2)
plot(px, pv)

plot(pc)

@assert Float32.(xcg2) ≈ Float32.(xcg1)

using Optim: ConjugateGradient, optimize

# recycle stuff from `state`
function g!(grad, x)
    for (Bxj, Bj) in zip(state.Bx, state.B)
        mul!(Bxj, Bj, x)
    end
    state.grad!(grad, state.npgrad, state.Bx)
    return grad
end

tmp = similar(x0)
@assert g!(state.grad_new, tmp) ≈ ∇f(tmp)

opt = optimize(f, g!, x0, ConjugateGradient())
xopt = opt.minimizer
@assert Float32.(xopt) ≈ Float32.(xcg1)

gui(); throw() # xx


#=
Now explore a fancier version
that needs less heap memory.
=#

work = LineSearchMMWork(uu, vv, α0) # pre-allocate
function lsmm1(gradf, curvf)
    return line_search_mm(gradf, curvf, uu, vv;
        ninner, out, fun, α0, work)
end
function lsmm2(dot_gradf, dot_curvf)
    gradn = [() -> nothing, () -> nothing]
    return line_search_mm(uu, vv, dot_gradf, dot_curvf;
        ninner, out, fun, α0, work)
end;

#=
The `let` statements below are a performance trick from the
[Julia manual](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured-1).
Using `Iterators.map` avoids allocating arrays like `z - y`
and does not even require any work space.
=#
gradz = [
    let y=y; z -> Iterators.map(-, z, y); end, # z - y
    let β=β, dpot=dpot; z -> Iterators.map(z -> β * dpot(z), z); end, # β * dψ.(z)
]
curvz = [
    1,
    let β=β, wpot=wpot; z -> Iterators.map(z -> β * wpot(z), z); end, # β * ωψ.(z)
]

function make_grad1c()
    w = similar(uu[1]) # work-space
    let w=w, y=y
        function grad1c(z)
            @. w = z - y
            return w
        end
    end
end

function make_grad2c()
    w = similar(uu[2]) # work-space
    let w=w, β=β, dpot=dpot
        function grad2c(z)
            @. w = β * dpot(z)
            return w
        end
    end
end

function make_curv2c()
    w = similar(uu[2]) # work-space
    let w=w, β=β, wpot=wpot
        function curv2c(z)
            @. w = β * wpot(z) # β * ωψ.(z)
            return w
        end
    end
end

gradc = [ # capture version
    make_grad1c(), # z - y
    make_grad2c(), # β * dψ.(z)
]
curvc = [
    1,
    make_curv2c(), # β * ωψ.(z)
]
#src @code_warntype curvc[2](vv[2]) # stable
#src @code_warntype gradc[1](vv[1]) # stable

sum_map(f::Function, args...) = sum(Iterators.map(f, args...))
dot_gradz = [
    let y=y; (v,z) -> sum_map((v,z,y) -> dot(v, z - y), v, z, y); end, # v'(z - y)
    let β=β, dpot=dpot; (v,z) -> β * sum_map((v,z) -> dot(v, dpot(z)), v, z); end, # β * (v'dψ.(z))
]
dot_curvz = [
    (v,z) -> norm(v)^2,
    let β=β, wpot=wpot; (v,z) -> β * sum_map((v,z) -> abs2(v) * wpot(z), v, z); end, # β * (abs2.(v)'ωψ.(z))
]
#src @code_warntype dot_gradz[1](vv[1], uu[1]) # stable with let
#src @code_warntype dot_gradz[2](vv[2], uu[2]) # stable with let
#src @code_warntype dot_curvz[2](vv[2], uu[2]) # stable with let

#src @btime dot_gradz[1]($(vv[1]), $(uu[1])) # 7 μs (1 allocation: 16 bytes)
#src @btime dot_gradz[2]($(vv[2]), $(uu[2])) # 1.9 μs (1 allocation: 16 bytes)
#src @btime dot_curvz[1]($(vv[1]), $(uu[1])) # 2. μs (1 allocation: 16 bytes)
#src @btime dot_curvz[2]($(vv[2]), $(uu[2])) # 1.9 μs (1 allocation: 16 bytes)

a1 = lsmm1(gradf, curvf)
a1c = lsmm1(gradc, curvc)
a2 = lsmm1(gradz, curvz)
a3 = lsmm2(dot_gradz, dot_curvz)
@assert a1 ≈ a2 ≈ a3 ≈ a1c

b1 = @benchmark a1 = lsmm1($gradf, $curvf)

#
bc = @benchmark a1c = lsmm1($gradc, $curvc)

#
b2 = @benchmark a2 = lsmm1($gradz, $curvz)

#
b3 = @benchmark a3 = lsmm2($dot_gradz, $dot_curvz)

#=
Timing results on my Mac:
- 95 μs
- 65 μs # 1c after using `make_`
- 80 μs
- 69 μs (and lowest memory)

The versions using `gradc` and `dot_gradz`
with their "properly captured" variables
are the fastest.
But all the versions here are pretty similar
so even using the simplest version
seems likely to be fine.
=#


#=
## Compare with LineSearches.jl

Was all this specialized effort useful?
Let's compare to the general line search methods in
[LineSearches.jl](https://github.com/JuliaNLSolvers/LineSearches.jl).

It seems that some of those methods do not allow ``α₀ = 0``
so we use 1.0 instead.
We use the default arguments for all the solvers,
which means some of them might terminate
before `ninner` iterations,
giving them a potential speed advantage.
=#

a0 = 1.0 # α0
hdh(α) = h(α), dh(α)
h0 = h(0)
dh0 = dh(0);
function ls_ls(linesearch)
    a1, fx = linesearch(h, dh, hdh, a0, h0, dh0)
    return a1
end;

solvers = [
    BackTracking( ; iterations = ninner),
    HagerZhang( ; linesearchmax = ninner),
    MoreThuente( ; maxfev = ninner),
]
for ls in solvers # check that they work properly
    als = ls_ls(ls)
    @assert isapprox(als, αstar; atol=1e-3)
end;

#
bbt = @benchmark ls_ls($(solvers[1]))

#
bhz = @benchmark ls_ls($(solvers[2]))

#
bmt = @benchmark ls_ls($(solvers[3]))

#=
On my Mac the timings are all much longer
compared to `line_search_mm`:
- 840 μs `BackTracking`
- 2.6 ms `HagerZhang`
- 3.9 ms `MoreThuente`

This comparison illustrates
the benefit of the "special purpose" line search.


The fastest version seems to be `BackTracking`,
so plot its iterates:
=#

alpha_bt = zeros(ninner + 1)
alpha_bt[1] = a0
for iter in 1:ninner
    tmp = BackTracking( ; iterations = iter)
    alpha_bt[iter+1] = ls_ls(tmp)
end
plot(0:ninner, alpha_bt, marker=:square, color=:blue,
    xlabel="Iteration", ylabel="BackTracking α")
plot!([0, ninner], [1,1] * αstar, color=:red)

#=
Unexpectedly,
`BackTracking` seems to terminate at the first iteration.
But even just that single iteration is slower than 7 iterations
of `line_search_mm`.
=#


#
prompt()


# ### Reproducibility

# This page was generated with the following version of Julia:

io = IOBuffer(); versioninfo(io); split(String(take!(io)), '\n')


# And with the following package versions

import Pkg; Pkg.status()
