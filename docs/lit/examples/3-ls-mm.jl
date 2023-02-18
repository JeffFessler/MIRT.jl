#=
# [Line search MM](@id 3-ls-mm)

Examples illustrating the
line-search method
based on majorize-minimize (MM) principles
in the Julia package
[`MIRT`](https://github.com/JeffFessler/MIRT.jl).
This method is probably most useful
for algorithm developers.
=#

#srcURL


# ### Setup

# Packages needed here.

using Plots; default(markerstrokecolor = :auto, label="")
using MIRTjim: prompt
using MIRT: line_search_mm, LineSearchMMWork
using LineSearches: BackTracking, HagerZhang, MoreThuente
using LinearAlgebra: norm, dot
using Random: seed!; seed!(0)
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
where ``D_j(z)`` is a positive semidefinite matrix
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
``f_1(u) = \frac{1}{2} \| u - y \|_2^2,``
``f_2(u) = β 1' ψ.(u),``
for which
``∇f_1(u) = u - y,``
``∇f_2(u) = β ψ.(u),``
and
``∇^2 f_1(u) = I,``
``∇^2 f_2(u) \succeq β \, \text{diag}(ω_{ψ}(u)).``

Set up an example and plot ``h(α)``.
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
fun(state) = state.α # log this
ninner = 7
out = Vector{Any}(undef, ninner+1)
α0 = 0
αstar = line_search_mm(gradf, curvf, uu, vv; ninner, out, fun, α0)
hmin = h(αstar)
scatter!([αstar], [hmin], marker=:star, color=:red)
scatter!([α0], [h(α0)], marker=:circle, color=:green)
ps = plot(0:ninner, out, marker=:circle, xlabel="iteration", ylabel="α",
    color = :green)
pd = plot(0:ninner, abs.(dh.(out)), marker=:diamond,
    yaxis = :log, color=:red,
    xlabel="iteration", ylabel="|dh(α)|")
pu = plot(1:ninner, log10.(max.(abs.(diff(out)), 1e-16)), marker=:square,
    color=:blue, xlabel="iteration", ylabel="log10(|α_k - α_{k-1}|)")
plot(pa, ps, pd, pu)

#=
Thanks to Huber's curvatures,
the ``α_t`` sequence converges very quickly.

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


include("../../../inc/reproduce.jl")
