#---------------------------------------------------------
# # [NUFFT](@id 2-nufft)
#---------------------------------------------------------

#=
Examples illustrating the `nufft` methods
in the Julia package
[`MIRT`](https://github.com/JeffFessler/MIRT.jl).

The `nufft` functions in this package
are wrappers around
[NFFT.jl](https://github.com/JuliaMath/NFFT.jl).

The plots on this page
are simply sanity checks
about of the approximation error
for that package
and they help verify the correctness
of the wrapper.

This page was generated from a single Julia file:
[2-nufft.jl](@__REPO_ROOT_URL__/2-nufft.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`2-nufft.ipynb`](@__NBVIEWER_ROOT_URL__/2-nufft.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`2-nufft.ipynb`](@__BINDER_ROOT_URL__/2-nufft.ipynb).


# ### Setup

# Packages needed here.

using Plots; default(markerstrokecolor = :auto, label="")
using MIRTjim: prompt
using MIRT: nufft_init, nufft_errors, dtft
using InteractiveUtils: versioninfo

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() && prompt(:prompt);


#=
## 1D example

This code illustrates how to call the NUFFT.
=#
Ω = range(-1,1,301) * π # digital frequencies (radians/sample)
N = 16 # # of samples
(nufft, nufft_adjoint, Anufft) = nufft_init(Ω, N)
x = randn(ComplexF32, N) # random 1D signal
y1 = nufft(x) # one way to compute NUFFT
y2 = Anufft * x # another way to compute NUFFT
@assert y1 == y2
yd = dtft(Ω, x) # exact (slow) nonuniform discrete Fourier transform
maximum(abs, yd - y1) # worst error for this 1D signal

# This plot shows that the NUFFT is very close to the exact nonuniform DFT.
xticks = ([-π,0,π], ["-π", "0", "π"]); xlims = (-π,π); xlabel="Ω"
p1 = plot(ylabel="Real part"; xlabel, xlims, xticks)
plot!(Ω, real(y1), label="real(NUFFT)", line=:dash)
plot!(Ω, real(yd), label="real(DTFT)", line=:dot)
p2 = plot(ylabel="Imag part"; xlabel, xlims, xticks)
plot!(Ω, imag(y1), label="imag(NUFFT)", line=:dash)
plot!(Ω, imag(yd), label="imag(NUFFT)", line=:dot)
p3 = plot(Ω, real(yd - y1), label="real(Error)"; xlabel, xlims, xticks)
p4 = plot(Ω, imag(yd - y1), label="imag(Error)"; xlabel, xlims, xticks)
plot(p1, p2, p3, p4)

#
prompt()

#=
## Plot worst-case errors vs ``N``

Plot worst-case error
over all frequencies ``Ω`` between ``0`` and ``2π/N``
for various ``N``.

Even for ``N=512``
the error is below `4e-6`.
=#
Nlist = 2 .^ (4:9)
Nlist = [Nlist..., 191, 385] # test odd N values
errfunn = N -> maximum(nufft_errors( ; N)[2])
elist = errfunn.(Nlist)
pn = scatter(Nlist, elist,
     xtick=Nlist, color=:red, xlabel="N", ylabel="error", yaxis=:log10)

#
prompt()


#=
## Plot error vs NFFT `m`.

For the default `m = 4`
the worst-case error is below `4e-6`.
=#
mlist = 3:7
errfunm = nfft_m -> maximum(nufft_errors( ; nfft_m)[2])
worst_m = errfunm.(mlist)
pm = scatter(mlist, worst_m,
    xlabel="m", ylabel="error", color=:magenta, yaxis=:log10)

#
prompt()


#=
## Plot error vs NFFT `sigma`.
For the default `σ = 4`
the worst-case error is below `4e-6`.
=#
slist = [1.5; 2:6] # σ
errfuns = nfft_sigma -> maximum(nufft_errors( ; nfft_sigma)[2])
worst_s = errfuns.(slist)
ps = scatter(slist, worst_s,
    xlabel="σ", ylabel="error", color=:blue, yaxis=:log10)

#
prompt()


# ### Reproducibility

# This page was generated with the following version of Julia:

io = IOBuffer(); versioninfo(io); split(String(take!(io)), '\n')


# And with the following package versions

import Pkg; Pkg.status()
