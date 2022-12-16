#---------------------------------------------------------
# # [Line search MM](@id 3-ls-mm)
#---------------------------------------------------------

#=
Examples illustrating the
line-search methods
based on majorize-minimize (MM) principles
in the Julia package
[`MIRT`](https://github.com/JeffFessler/MIRT.jl).

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
using MIRT: ls_mm
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
where ``\mathbb{F}`` is the field of real or complex numbers,
matrix ``B_j`` has size ``M_j × N``,
and ``f_j : \mathbb{F}^{M_j} ↦ \mathbb{R}``.

Of course one could try to apply general-purpose optimization methods,
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
α_* = \arg\min_{α in \mathbb{R}} h(α)
,\qquad
h(α) = f(x + α d),
```
for some search direction `d`.

There are general purpose line search algorithms in
[LineSearches.jl](https://github.com/JuliaNLSolvers/LineSearches.jl)
but here we focus on the specific form given above.

For that form we see that we have the special structure
```math
h(α) = \sum_{j=1}^J f_j(u_j + α v_j)
,\qquad
u_j = B_j x
,\quad
v_j = B_j d
.```

Here we focus further
on the case where each function ``f_j(⋅)``
has a quadratic majorizer
of the form
```math
f_j(x) ≤ q_j(x,z) = f_j(z) + \text{real}(⟨ ∇f_j(z), x - z ⟩)
+ \frac{1}{2} (x - z)' D_j(z) (x - z)
,```
where typically ``D_j(z)`` is a diagonal matrix.
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
α_{t+1} = α_t - c_1(α_t) / c_2(α_t).
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
=#

#
prompt()


# ### Reproducibility

# This page was generated with the following version of Julia:

io = IOBuffer(); versioninfo(io); split(String(take!(io)), '\n')


# And with the following package versions

import Pkg; Pkg.status()
