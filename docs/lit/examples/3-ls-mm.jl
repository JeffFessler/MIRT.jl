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
using MIRT: line_search_mm, LineSearchMMWork
using LinearAlgebra: norm, dot
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

using BenchmarkTools: @btime # todo ndgrid

# Fair potential, its derivative and Huber weighting function
function fair_pot()
    fpot(z,δ) = δ^2 * (abs(z/δ) - log(1 + abs(z/δ)))
    dpot(z,δ) = z / (1 + abs(z/δ))
    wpot(z,δ) = 1 / (1 + abs(z/δ))
    return fpot, dpot, wpot
end

M, N = 8000, 1000
A = randn(M,N)
x0 = randn(N) .* (rand(N) .< 0.4) # sparse vector
y = A * x0 + 0.001 * randn(M)
β = 5
δ = 0.1
#ψ, dψ, ωψ = ((d) -> fair_pot(d))(δ) # dirty trick!?
#ψ, dψ, ωψ = Base.Fix2.(fair_pot(), δ)
fpot, dpot, wpot = Base.Fix2.(fair_pot(), δ)
#@code_warntype wpot(5.)
#throw()

#ψ(z) = δ^2 * (abs(z/δ) - log(1 + abs(z/δ)))
#dψ(z) = z / (1 + abs(z/δ))
#ωψ(z) = 1 / (1 + abs(z/δ))
#ψ1 = (z) -> ψ2(z, δ)
#dψ1 = (z) -> dψ2(z, δ)
#ωψ1 = (z) -> ωψ2(z, δ)

#ψ2(z,δ) = δ^2 * (abs(z/δ) - log(1 + abs(z/δ)))
dψ2(z,δ) = z / (1 + abs(z/δ))
ωψ2(z,δ) = 1 / (1 + abs(z/δ))

f(x) = 0.5 * norm(A * x - y)^2 + β * sum(fpot, x)
∇f(x) = A' * (A * x - y) + β * dpot.(x)
x = randn(N) # random point
d = -∇f(x)/40 # some search direction
h(α) = f(x + α * d)
dh(α) = d' * ∇f(x + α * d)
pa = plot(h, xlabel="α", ylabel="h(α)", xlims=(-1, 2))


# Apply MM-based line search: simple version
gradf = [ # key ingredients!
    u -> u - y,
    u -> β * dpot.(u), # slower
#   u -> β * dψ2.(u, δ), # faster!?
]
curvf = [1, u -> β * wpot.(u)] # slower
#curvf = [1, u -> β * ωψ2.(u,δ)] # faster!

uu = [A * x, x]
vv = [A * d, d]
fun(state, iter) = state.α
ninner = 7
out = Vector{Any}(undef, ninner+1)
out .= 0 # todo
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

work = LineSearchMMWork(uu, vv, α0) # pre-allocate
function lsmm(gradf, curvf)
    return line_search_mm(uu, vv, gradf, curvf;
        ninner, out, fun, α0, work)
end
gradn = [() -> nothing, () -> nothing]
function lsmm2(dot_gradf, dot_curvf)
    return line_search_mm(uu, vv, gradn, gradn;
        ninner, out, fun, α0, work, dot_gradf, dot_curvf)
end

precompile(dpot, (Float64,))

# faster? version of key ingredients! todo time
gradz = [
    z -> Iterators.map(-, z, y), # z - y
    z -> Iterators.map((z, β, δ) -> β * dψ2(z, δ), z, # β * dψ.(z)
        Iterators.cycle(β),
        Iterators.cycle(δ))
#   z -> Iterators.map((z, β) -> β * dpot(z), z, # β * dψ.(z)
#       Iterators.cycle(β))
]
curvz = [
    1,
#   z -> Iterators.map((z, β, δ) -> β * ωψ2(z, δ), z, # β * ωψ.(z)
#       Iterators.cycle(β),
#       Iterators.cycle(δ))
    z -> Iterators.map((z, β) -> β * wpot(z), z, # β * ωψ.(z)
        Iterators.cycle(β))
]

sum_map(f::Function, args...) = sum(Iterators.map(f, args...))
dot_gradz = [
    (v,z) -> sum_map((v,z,y) -> dot(v, z - y), v, z, y), # v'(z - y)
    (v,z) -> β * sum_map((v,z,d) -> dot(v, dψ2(z, d)), v, z,
        Iterators.cycle(δ)), # β * (v'dψ.(z))
#   (v,z) -> β * sum_map((v,z) -> dot(v, dpot(z)), v, z), # β * (v'dψ.(z))
]
dot_curvz = [
    (v,z) -> norm(v)^2,
    (v,z) -> β * sum_map((v,z,d) -> abs2(v) * ωψ2(z, d), v, z,
        Iterators.cycle(δ)), # β * (abs2.(v)'ωψ.(z))
#   (v,z) -> β * sum_map((v,z) -> abs2(v) * wpot(z), v, z), # β * (abs2.(v)'ωψ.(z))
]

# @btime dψ(4,$δ) # 2.928 ns (0 allocations: 0 bytes)
# @btime dψ(4) # 66.980 ns (4 allocations: 64 bytes) todo why?  closure?
#tmp = dot_gradz2(vv[2], uu[2]) # warm-up
@btime dot_gradz[1]($(vv[1]), $(uu[1])) # 2.6 μs (7 allocations: 160 bytes)
@btime dot_gradz[2]($(vv[2]), $(uu[2])) # 2.6 μs (7 allocations: 160 bytes)

@btime dot_curvz[1]($(vv[1]), $(uu[1]))
@btime dot_curvz[2]($(vv[2]), $(uu[2])) # 2.8 μs (9 allocations: 208 bytes)
throw()

a1 = lsmm(gradf, curvf)
a2 = lsmm(gradz, curvz)
a3 = lsmm2(dot_gradz, dot_curvz)
@assert a1 ≈ a2 ≈ a3
#throw()

# todo timing comparisons, fancier use
@btime a1 = lsmm($gradf, $curvf)
@btime a2 = lsmm($gradz, $curvz)
@btime a3 = lsmm2($dot_gradz, $dot_curvz)

#αplot =
#
# prompt()

# ### Reproducibility

# This page was generated with the following version of Julia:

# io = IOBuffer(); versioninfo(io); split(String(take!(io)), '\n')


# And with the following package versions

# import Pkg; Pkg.status() # todo


#=
#junk

#d1(v,z,y) = v * (z - y) # need conj in complex case
#d2(v,z) = v * dψ(z)
#d2 = (v,z) -> v * dψ1(z)
#d2(v,z) = v * dψ(z,0.1) # fast
#d2(v,z) = v * dψ(z,copy(δ))
#d2(v,z) = v * dψ(z,eval(:δ)) # slow
#d2(v,z,δ) = v * dψ(z,δ) # todo try?
#make_d2(d) = (v,z) -> v * dψ(z, d) # fast
#maker_d2 = (d) -> (v,z) -> v * dψ(z, d) # fast
#maker_d2 = (v,z) -> v * dψ(z) # slow
#maker_d2 = () -> (v,z) -> v * dψ(z) # slow
#maker_d2 = () -> (v,z) -> v * dψ(z,δ) # slow
#d2 = maker_d2()
#d2 = ((d) -> (v,z) -> v * dψ2(z, d))(δ) # fast - a dirty trick!?
#d2 = (() -> (v,z) -> v * dψ(z))() # slow

function mydot1(v, z) # slow with big allocations - why?
     total = zero(v[1] * (z[1] - y[1]))
     for i in eachindex(v)
         total += v[i] * (z[i] - y[i])
     end
     return total
end

function mydot2(v, z)
     total = zero(v[1] * dψ(z[1]))
     for i in eachindex(v)
         total += v[i] * dψ(z[i])
     end
     return β * total
end
=#

#=

function fair_pot(δ)
    ψ = (z) -> δ^2 * (abs(z/δ) - log(1 + abs(z/δ)))
    dψ = (z) -> z / (1 + abs(z/δ))
    ωψ = (z) -> 1 / (1 + abs(z/δ))

    ψ = (z,d) -> d^2 * (abs(z/d) - log(1 + abs(z/d)))
    dψ = (z,d) -> z / (1 + abs(z/d))
    ωψ = (z,d) -> 1 / (1 + abs(z/d))

    dψ = ((d) -> ((z) -> dψ(z,d)))(δ)

    ψ(z) = δ^2 * (abs(z/δ) - log(1 + abs(z/δ)))
    dψ(z) = z / (1 + abs(z/δ))
    ωψ(z) = 1 / (1 + abs(z/δ))

    return ψ, dψ, ωψ
end

ψ, dψ, ωψ = fair_pot(δ)
ψ, dψ, ωψ = fair_pot(0.1)
ψ2, dψ2old, ωψ2old = fair_pot()
=#
