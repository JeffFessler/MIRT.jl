# Standalone experiment, deliberately not registered in the CI benchmark suite.
# Run from any directory: julia --startup-file=no /path/to/MIRT.jl/benchmark/pogm_restart_gr.jl
# Optional arguments: samples (default 5), rows (600), columns (60).
# Loads solver SOURCE from this checkout; does not use an installed MIRT package.
# Uses only Julia standard libraries, so no package installation is necessary.
module POGMRestartBenchmark

using LinearAlgebra: BLAS, norm
using Statistics: median
using Printf: @printf

include(joinpath(@__DIR__, "../src/algorithm/general/pogm_restart.jl"))
include(joinpath(@__DIR__, "../src/algorithm/general/pogm_restart_gr.jl"))
include(joinpath(@__DIR__, "../test/algorithm/general/pogm_test_problems.jl"))

function solve(p, mode, budget, tolerance)
    options = (; niter=budget, g_prox=p.prox)
    if mode == :original
        return first(pogm_restart(p.x0, p.cost, p.grad, p.L; options...))
    elseif mode in (:bang_fixed, :bang_check5)
        stopping = mode == :bang_check5 ? (; atol=tolerance, check_every=5) : (;)
        return first(pogm_restart_gr(p.x0, p.cost, p.grad, p.L; options...,
            f_grad! = p.grad!, g_prox! = p.prox!, stopping...))
    elseif mode == :fixed
        return first(pogm_restart_gr(p.x0, p.cost, p.grad, p.L; options...))
    else
        return first(pogm_restart_gr(p.x0, p.cost, p.grad, p.L;
            options..., atol=tolerance, check_every=mode))
    end
end

# Counts are collected separately, outside the timed runs and accuracy checks.
function diagnostics(p, mode, budget, tolerance)
    gradients, proximals, costs = Ref(0), Ref(0), Ref(0)
    grad = x -> (gradients[] += 1; p.grad(x))
    prox = (z, c) -> (proximals[] += 1; p.prox(z, c))
    cost = x -> (costs[] += 1; p.cost(x))
    if mode == :original
        x, _ = pogm_restart(p.x0, cost, grad, p.L; niter=budget, g_prox=prox)
        iterations, reason = budget, :niter
    else
        options = if mode in (:bang_fixed, :bang_check5)
            grad! = (d,x) -> (gradients[] += 1; p.grad!(d,x))
            prox! = (d,z,c) -> (proximals[] += 1; p.prox!(d,z,c))
            stopping = mode == :bang_check5 ? (; atol=tolerance, check_every=5) : (;)
            (; f_grad! = grad!, g_prox! = prox!, stopping...)
        else
            mode == :fixed ? (;) : (; atol=tolerance, check_every=mode)
        end
        x, _, info = pogm_restart_gr(p.x0, cost, grad, p.L;
            niter=budget, g_prox=prox, return_info=true, options...)
        iterations, reason = info.iterations, info.reason
    end
    return (; x, iterations, reason, gradients=gradients[],
        proximals=proximals[], costs=costs[])
end

# Cheap gradient (one result vector), identity prox, allocation-free scalar cost:
# this exposes solver temporaries rather than dense loss-function work.
function allocation_probe()
    x = ones(1024)
    grad = x -> 0.25*x
    cost = x -> sum(abs2, x)/8
    println("Allocation probe: Float64[1024], one-vector gradient, identity prox, no stopping.")
    for (label, solver) in (("original", pogm_restart), ("isolated", pogm_restart_gr))
        totals = Int[]
        for budget in (100, 200)
            solver(x, cost, grad, 1.; niter=budget) # warm
            bytes = [@allocated(solver(x, cost, grad, 1.; niter=budget)) for _ in 1:3]
            push!(totals, Int(median(bytes)))
        end
        println(label, ": bytes(100)=", totals[1], "; bytes(200)=", totals[2],
            "; incremental bytes/iteration=", (totals[2]-totals[1])/100)
    end
    println("Increment includes gradient results and history growth; not a claim of zero hot-loop allocations.\n")
end

function callback_allocations(p)
    dest = similar(p.x0)
    p.grad(p.x0); p.grad!(dest, p.x0)
    p.prox(p.x0, 0.3); p.prox!(dest, p.x0, 0.3)
    return (; gradient=@allocated(p.grad(p.x0)),
        gradient_bang=@allocated(p.grad!(dest, p.x0)),
        prox=@allocated(p.prox(p.x0, 0.3)),
        prox_bang=@allocated(p.prox!(dest, p.x0, 0.3)))
end

function main(; samples=5, m=600, n=60)
    samples >= 1 || throw(ArgumentError("samples must be positive"))
    m >= n >= 5 || throw(ArgumentError("require rows >= columns >= 5"))
    # Avoid small matrix timings dominated by multithreaded BLAS overhead.
    previous_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        println("Julia ", VERSION, "; BLAS threads=1; samples=", samples,
            "; dimensions=", (m, n))
        println("Source: ", realpath(joinpath(@__DIR__, "../src/algorithm/general/pogm_restart_gr.jl")))
        println("Seed=17; lambda=0.015; budget=2000; absolute residual tolerance=1e-7")
        println("Median elapsed time and allocated bytes over warmed full calls, including initialization/history.")
        println("GC is run before each sample (outside timing); GC during a solve is included.")
        println("Counts exclude setup/reference/diagnostics. No callback snapshots are retained.")
        println("This measures complete calls, not hot-loop-only allocations. No speed assertions.\n")
        allocation_probe()
        budget, tolerance = 2000, 1e-7
        for kind in (:logistic, :squared_hinge)
            p = pogm_test_problem(kind; m, n)
            println(kind, " callback allocated bytes (scratch initialization excluded): ",
                callback_allocations(p))
            reference, _ = pogm_restart(p.x0, p.cost, p.grad, p.L;
                g_prox=p.prox, niter=10000)
            reference_residual = p.residual(reference)
            reference_residual < 1e-10 || error("reference did not reach required accuracy")
            reference_cost = p.cost(reference)
            println(kind, "; reference residual=", reference_residual,
                "; reference cost=", reference_cost)
            # Offline oracle: find a fixed budget meeting the SAME accuracy gates.
            # Calibration is excluded from timings, and is not a deployable stopping rule.
            calibrated = nothing
            for k in 1:budget
                candidate = solve(p, :fixed, k, tolerance)
                if p.residual(candidate) <= tolerance &&
                    abs(p.cost(candidate)-reference_cost) < 1e-8 &&
                    norm(candidate-reference)/max(1., norm(reference)) < 1e-4
                    calibrated = k
                    break
                end
            end
            calibrated === nothing && error("no fixed budget reached target accuracy")
            println("Offline calibrated fixed budget=", calibrated,
                " (calibration excluded; requires advance knowledge).")
            println("mode          ms       MiB    iterations   grad   prox   cost    residual      cost_gap      rel_x_error   accuracy_ok  reason")
            cases = ((:original, budget, "original"), (:fixed, budget, "fixed"),
                (:original, calibrated, "orig_match"), (:fixed, calibrated, "fixed_match"),
                (1, budget, "check_1"), (5, budget, "check_5"), (20, budget, "check_20"),
                (:bang_fixed, budget, "bang_fixed"), (:bang_check5, budget, "bang_check5"))
            for (mode, run_budget, label) in cases
                solve(p, mode, run_budget, tolerance) # compile and warm this path
                elapsed, allocated = Float64[], Int[]
                for _ in 1:samples
                    GC.gc()
                    sample = @timed solve(p, mode, run_budget, tolerance)
                    push!(elapsed, sample.time)
                    push!(allocated, sample.bytes)
                end
                info = diagnostics(p, mode, run_budget, tolerance)
                residual = p.residual(info.x)
                gap = p.cost(info.x) - reference_cost
                relative_error = norm(info.x-reference) / max(1., norm(reference))
                accuracy_ok = residual <= tolerance && abs(gap) < 1e-8 && relative_error < 1e-4
                @printf("%-11s %8.3f %9.3f %8d %7d %6d %6d %12.3e %12.3e %12.3e %8s  %s\n",
                    label, 1000*median(elapsed), median(allocated)/2.0^20,
                    info.iterations, info.gradients, info.proximals, info.costs,
                    residual, gap, relative_error, string(accuracy_ok), string(info.reason))
            end
            println()
        end
    finally
        BLAS.set_num_threads(previous_threads)
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    values = parse.(Int, ARGS)
    length(values) <= 3 || error("usage: pogm_restart_gr.jl [samples [rows [columns]]]")
    POGMRestartBenchmark.main(;
        samples=length(values) >= 1 ? values[1] : 5,
        m=length(values) >= 2 ? values[2] : 600,
        n=length(values) >= 3 ? values[3] : 60)
end
