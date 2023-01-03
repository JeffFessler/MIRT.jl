# test/runtests.jl

import MIRT
using Plots: default, plot
using Test: @testset, detect_ambiguities

if !isinteractive() # for usual noninteractive test, do not show plots
    default(show=false, reuse=true) # trick from Plots/runtests.jl
    ENV["GKSwstype"] = "100" # also from Plots/runtests.jl
end


"""
    mirt_setup_plot()

MIRT tests use Plots (with default GR backend).
When running on a remote linux server over X11 I get error messages from GKS.
The following kludge (?) tries to overcome this problem.
"""
function mirt_setup_plot()

    try
        plot(ones(3))
    catch
        if Sys.islinux
            ENV["GKS_WSTYPE"] = "x11"
        else
            @warn "Plot not working for Sys.KERNEL=$(Sys.KERNEL); test may fail"
        end
    end

    try
        plot(ones(3))
    catch
        @warn "Plot still not working; test may fail"
    end

    nothing
end

mirt_setup_plot()


# Run all MIRT tests.
# To avoid plotting during tests, first do `default(show=false)`

otherlist = (
"../data",
)

basedir = dirname(pathof(MIRT))
for root in otherlist
    @testset "$root" begin
        tmp = joinpath(basedir, root, "z-test.jl")
        include(tmp)
    end
end


srclist = [
"algorithm",
"fbp",
"io",
"mri",
"nufft",
"regularize",
"system",
"utility",
]

for root in srclist
    tmp = joinpath(root, "z-test.jl")
#   printstyled(tmp * "\n", color=:blue, bold=true)
    include(tmp)
end

@test length(detect_ambiguities(MIRT)) == 0
