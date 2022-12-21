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

M, N = 1000, 2000
A = randn(M,N)
x0 = randn(N) .* (rand(N) .< 0.4) # sparse vector
y = A * x0 + 0.001 * randn(M)
β = 95
δ = 0.1
fpot, dpot, wpot = Base.Fix2.(fair_pot(), δ)

f(x) = 0.5 * norm(A * x - y)^2 + β * sum(fpot, x)
∇f(x) = A' * (A * x - y) + β * dpot.(x)
x = randn(N) # random point
d = -∇f(x)/M # some search direction
h(α) = f(x + α * d)
dh(α) = d' * ∇f(x + α * d)
pa = plot(h, xlabel="α", ylabel="h(α)", xlims=(-1, 2))


# Apply MM-based line search: simple version.
# The key inputs are the gradient and curvature functions:
gradf = [
    u -> u - y, # ∇f₁ for data-fit term
    u -> β * dpot.(u), # ∇f₂ for regularizer
]
curvf = [
    1, # curvature for data-fit term
    u -> β * wpot.(u), # Huber curvature for regularizer
]

uu = [A * x, x] # [u₁ u₂]
vv = [A * d, d] # [v₁ v₂]
fun(state, iter) = state.α # log this
ninner = 7
out = Vector{Any}(undef, ninner+1)
#out .= 0 # todo
α0 = 0
αstar = line_search_mm(uu, vv, gradf, curvf; ninner, out, fun, α0)
scatter!([αstar], [h(αstar)], marker=:star, color=:red)
scatter!([α0], [h(α0)], marker=:circle, color=:green)
ps = plot(0:ninner, out, marker=:circle, xlabel="iteration", ylabel="α",
    color = :green)
pd = plot(0:ninner, abs.(dh.(out)), marker=:diamond,
    yaxis = :log, color=:red,
    xlabel="iteration", ylabel="|dh(α)|")
pu = plot(1:ninner, log10.(max.(abs.(diff(out)), 1e-16)), marker=:square,
    color=:blue, xlabel="iteration", ylabel="log10(|α_k - α_{k-1}|)")
plot(pa, ps, pd, pu)
# gui(); throw()

# Thanks to Huber's curvatures,
# the ``α`` sequence converges very quickly.

work = LineSearchMMWork(uu, vv, α0) # pre-allocate
function lsmm1(gradf, curvf)
    return line_search_mm(uu, vv, gradf, curvf;
        ninner, out, fun, α0, work)
end
function lsmm2(dot_gradf, dot_curvf)
    gradn = [() -> nothing, () -> nothing]
    return line_search_mm(uu, vv, gradn, gradn;
        ninner, out, fun, α0, work, dot_gradf, dot_curvf)
end

# The `let` statements here are a performance trick from
# https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured-1
# Using `Iterators.map` avoids allocating arrays like `z - y`
# and does not even require any work space
gradz = [
    let y=y; z -> Iterators.map(-, z, y); end, # z - y
    let β=β, dpot=dpot; z -> Iterators.map(z -> β * dpot(z), z); end, # β * dψ.(z)
]
curvz = [
    1,
    let β=β, wpot=wpot; z -> Iterators.map(z -> β * wpot(z), z); end, # β * ωψ.(z)
]

#=
ww = similar.(uu)
function grad1c(z,y)
    w = ww[1]
    @. w = z - y
    return w
end

function grad2c(z,β,dpot)
    w = ww[2]
    @. w = β * dpot(z)
    return w
end

function curv2c(z,β,wpot)
    w = ww[2]
    @. w = β * wpot(z)
    return w
end
=#

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
#   let y=y; z -> grad1c(z, y); end, # z - y
    make_grad1c(), # z - y
    make_grad2c(), # β * dψ.(z)
#   let β=β, dpot=dpot; z -> grad2c(z, β, dpot); end, # β * dψ.(z)
]
curvc = [
    1,
#   let β=β, wpot=wpot; z -> curv2c(z, β, wpot); end, # β * ωψ.(z)
    make_curv2c(), # β * ωψ.(z)
]
#@code_warntype curvc[2](vv[2]) # stable
#@code_warntype gradc[1](vv[1]) # stable

sum_map(f::Function, args...) = sum(Iterators.map(f, args...))
dot_gradz = [
    let y=y; (v,z) -> sum_map((v,z,y) -> dot(v, z - y), v, z, y); end, # v'(z - y)
    let β=β, dpot=dpot; (v,z) -> β * sum_map((v,z) -> dot(v, dpot(z)), v, z); end, # β * (v'dψ.(z))
]
dot_curvz = [
    (v,z) -> norm(v)^2,
    let β=β, wpot=wpot; (v,z) -> β * sum_map((v,z) -> abs2(v) * wpot(z), v, z); end, # β * (abs2.(v)'ωψ.(z))
]
#@code_warntype dot_gradz[1](vv[1], uu[1]) # stable with let
#@code_warntype dot_gradz[2](vv[2], uu[2]) # stable with let
#@code_warntype dot_curvz[2](vv[2], uu[2]) # stable with let

#@btime dot_gradz[1]($(vv[1]), $(uu[1])) # 7 μs (1 allocation: 16 bytes)
#@btime dot_gradz[2]($(vv[2]), $(uu[2])) # 1.9 μs (1 allocation: 16 bytes)
#@btime dot_curvz[1]($(vv[1]), $(uu[1])) # 2. μs (1 allocation: 16 bytes)
#@btime dot_curvz[2]($(vv[2]), $(uu[2])) # 1.9 μs (1 allocation: 16 bytes)

a1 = lsmm1(gradf, curvf)
a1c = lsmm1(gradc, curvc)
a2 = lsmm1(gradz, curvz)
a3 = lsmm2(dot_gradz, dot_curvz)
@assert a1 ≈ a2 ≈ a3 ≈ a1c

@btime a1 = lsmm1($gradf, $curvf)
@btime a1c = lsmm1($gradc, $curvc)
@btime a2 = lsmm1($gradz, $curvz)
@btime a3 = lsmm2($dot_gradz, $dot_curvz)

#=
Timing results on my mac:
127.552 μs (382 allocations: 505.69 KiB)
117.703 μs (421 allocations: 11.66 KiB) # 1c before using using make_
 89.720 μs (319 allocations: 9.41 KiB) # 1c after using make_ !!
103.408 μs (316 allocations: 8.81 KiB)
92.245 μs (233 allocations: 5.06 KiB)
=#

# todo use "let" within make_dot_ constructors!

# todo: compare with LS package


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
