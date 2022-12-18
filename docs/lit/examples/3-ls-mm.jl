#=
# [Line search MM](@id 3-ls-mm)

Examples illustrating the
line-search method
based on majorize-minimize (MM) principles
in the Julia package
[`MIRT`](https://github.com/JeffFessler/MIRT.jl).
This method is probably most useful
for algorithm developers.

This page was generated from a single Julia file:
[3-ls-mm.jl](@__REPO_ROOT_URL__/3-ls-mm.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`3-ls-mm.ipynb`](@__NBVIEWER_ROOT_URL__/3-ls-mm.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`3-ls-mm.ipynb`](@__BINDER_ROOT_URL__/3-ls-mm.ipynb).


# ### Setup

# Packages needed here.

using Plots; default(markerstrokecolor = :auto, label="")
using MIRTjim: prompt
using MIRT: line_search_mm
using LinearAlgebra: norm
using Random: seed!; seed!(0)
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

Many algorithms for solving such problems
require an inner 1D optimization problem
called a
[line search](https://en.wikipedia.org/wiki/Line_search)
of the form
```math
α_* = \arg\min_{α ∈ \mathbb{R}} h(α)
,\qquad
h(α) = f(x + α d),
```
for some search direction ``d``.
There are general purpose line search algorithms in
[LineSearches.jl](https://github.com/JuliaNLSolvers/LineSearches.jl),
but here we focus on the specific form of ``f`` given above.

For that form we see that we have the special structure
```math
h(α) = \sum_{j=1}^J f_j(u_j + α v_j)
,\qquad
u_j = B_j x
,\quad
v_j = B_j d.
```

Here we focus further
on the case where each function ``f_j(⋅)``
has a quadratic majorizer
of the form
```math
f_j(x) ≤ q_j(x,z) = f_j(z) + \text{real}(⟨ ∇f_j(z), x - z ⟩)
+ \frac{1}{2} (x - z)' D_j(z) (x - z),
```
where ``D_j(z)`` is a positive definite matrix
that typically is diagonal.
Often it is a constant times the identity matrix,
e.g.,
a Lipschitz constant for ``∇f_j``,
but often there are
sharper majorizers.

Such quadratic majorizers
induce a quadratic majorizer
for ``h(α)`` as well:
```math
h(α) ≤ q(α; α_t) =
\sum_{j=1}^J q_j(u_j + α v_j; u_j + α_t v_j)
= h(α_t) + c_1(α_t) (α - α_t)
+ \frac{1}{2} c_2(α_t) (α - α_t)^2
```
where
```math
c_1(α_t) = \sum_{j=1}^J \text{real}(⟨ ∇f_j(u_j + α_t v_j), v_j ⟩)
,\qquad
c_2(α_t) = \sum_{j=1}^J v_j' D_j(u_j + α_t v_j) v_j.
```

The `line_search_mm` function
in this package
uses this quadratic majorizer
to update ``α``
using the iteration
```math
α_{t+1}
= \arg\min_{α} q(α;α_t)
= α_t - c_1(α_t) / c_2(α_t).
```
Being an MM algorithm,
it is guaranteed to decrease
``h(α)`` every update.
For an early exposition of this approach,
see
[Fessler & Booth, 1999](http://doi.org/10.1109/83.760336).

From the above derivation,
the main ingredients needed are
functions for computing
the dot products
``⟨ ∇f_j(u_j + α_t v_j), v_j ⟩``
and
``v_j' D_j(u_j + α_t v_j) v_j``.

The `line_search_mm` function
can construct such functions
given input gradient functions
``[∇f_1,…,∇f_J]``
and curvature functions
``[ω_1,…,ω_J]``
where
``D_j(z) = \text{Diag}(ω_j(z))``.

Alternatively,
the user can provide functions
for computing the dot products.

All of this is best illustrated
by an example.

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
``f_1(u) = \frac{1}{2} \| u - y \|_2^2``,
``f_2(u) = β 1' ψ.(u)``,
for which
``∇f_1(u) = u - y``,
``∇f_2(u) = β ψ.(u)``,
and
``∇^2 f_1(u) = I``,
``∇^2 f_2(u) \succeq β \text{diag}(ω_{ψ}(u))``.

Set up an example and plot ``h(α)``.
=#

M, N = 20, 10
A = randn(M,N)
x0 = randn(N) .* (rand(N) .< 0.4) # sparse vector
y = A * x0 + 0.001 * randn(M)
β = 5
δ = 0.1
ψ(z) = δ^2 * (abs(z/δ) - log(1 + abs(z/δ)))
dψ(z) = z / (1 + abs(z/δ))
ωψ(z) = 1 / (1 + abs(z/δ))
f(x) = 0.5 * norm(A * x - y)^2 + β * sum(ψ, x)
∇f(x) = A' * (A * x - y) + β * dψ.(x)
x = randn(N) # random point
d = -∇f(x)/40 # some search direction
h(α) = f(x + α * d)
dh(α) = d' * ∇f(x + α * d)
pa = plot(h, xlabel="α", ylabel="h(α)", xlims=(-1, 3), ylims=(0,100))


# Apply MM-based line search: simple version
gradf = [u -> u - y, u -> β * dψ.(u)] # key ingredients!
curvf = [1, u -> β * ωψ.(u)]

# faster? version of key ingredients! todo time
gradf = [
    u -> Iterators.map(Base.splat(-), zip(u,y)), # u - y
    u -> Iterators.map(Base.Fix1(*,β) ∘ dψ, u), # β * dψ.(u)
]
curvf = [
    1,
    u -> Iterators.map(Base.Fix1(*,β) ∘ ωψ, u), # β * ωψ.(u)
]

uu = [A * x, x]
vv = [A * d, d]
fun(state, iter) = state.α
ninner = 7
out = Vector{Any}(undef, ninner+1)
α0 = 0
αstar = line_search_mm(uu, vv, gradf, curvf; ninner, out, fun, α0)
scatter!([αstar], [h(αstar)], marker=:star, color=:red)
scatter!([α0], [h(α0)], marker=:circle, color=:green)
ps = plot(0:ninner, out, marker=:circle, xlabel="iteration", ylabel="α",
    color = :green)
pd = plot(0:ninner, abs.(dh.(out)), marker=:diamond,
    yaxis = :log, color=:red,
    xlabel="iteration", ylabel="|dh(α)|")
plot(pa, ps, pd)

# Thanks to Huber's curvatures,
# the ``α`` sequence converges very quickly.

# todo timing comparisons, fancier use

#αplot =
#
# prompt()

# ### Reproducibility

# This page was generated with the following version of Julia:

# io = IOBuffer(); versioninfo(io); split(String(take!(io)), '\n')


# And with the following package versions

# import Pkg; Pkg.status() # todo
