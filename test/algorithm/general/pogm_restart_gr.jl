# Compare the full trajectory with the original POGM + gradient restart path.
using MIRT: pogm_restart, pogm_restart_gr
using LinearAlgebra: norm
using Test: @test, @testset, @test_throws

@testset "POGM gradient restart isolation" begin
    # Diagonal quadratic: exact Lipschitz constant, deterministic and ill-conditioned.
    d = [1., 3., 20., 100.]
    b = [2., -1., 3., -2.]
    x0 = [4., -3., 2., 1.]
    f_grad = x -> d .* x - b
    smooth = x -> sum(d .* abs2.(x)) / 2 - sum(b .* x)
    snapshot = (iter, x, y, restarted) -> (iter, copy(x), copy(y), restarted)
    saw_restart = false
    saw_no_restart = false

    for penalty in (:smooth, :l1, :l2), mu in (0., 1.),
        bsig in (1., 0.5), cutoff in (0., 0.4), niter in (0, 1, 2, 80)
        prox = penalty == :l1 ? ((z, c) -> sign.(z) .* max.(abs.(z) .- 0.3*c, 0)) :
            penalty == :l2 ? ((z, c) -> z / (1 + 0.3*c)) : ((z, c) -> z)
        cost = x -> smooth(x) + (penalty == :l1 ? 0.3*sum(abs, x) :
            penalty == :l2 ? 0.3/2*sum(abs2, x) : 0.)
        old_calls = Ref(0)
        new_calls = Ref(0)
        old_cost = x -> (old_calls[] += 1; cost(x))
        new_cost = x -> (new_calls[] += 1; cost(x))
        options = (; f_mu=mu, bsig, restart_cutoff=cutoff, niter,
            g_prox=prox, fun=snapshot)
        xref, href = pogm_restart(x0, old_cost, f_grad, 100.;
            mom=:pogm, restart=:gr, options...)
        x, history = pogm_restart_gr(x0, new_cost, f_grad, 100.; options...)
        @test x == xref
        @test history == href # both iterates and restart flags, every iteration
        @test length(history) == niter + 1
        @test old_calls[] == new_calls[] == niter + 1
        if niter > 0
            @test x == history[end][2] # secondary iterate is returned
            saw_restart |= any(h -> h[4], history[2:end])
            saw_no_restart |= any(h -> !h[4], history[2:end])
        end
    end
    @test saw_restart
    @test saw_no_restart

    # Exercise default prox and callback, and inherited argument checks.
    x, out = pogm_restart_gr(x0, smooth, f_grad, 100.; niter=1)
    xref, oref = pogm_restart(x0, smooth, f_grad, 100.; niter=1)
    @test x == xref
    @test out == oref
    for options in ((; f_mu=-1.), (; bsig=-1.), (; restart_cutoff=1.),
        (; restart_cutoff=-1.))
        @test_throws ArgumentError pogm_restart_gr(x0, smooth, f_grad, 100.; options...)
    end
    @test_throws ArgumentError pogm_restart_gr(x0, smooth, f_grad, -1.)
end

@testset "POGM early stopping" begin
    d = [1., 3., 20., 100.]
    b = [0.1, -1., 3., -2.]
    x0 = [4., -3., 2., 1.]
    lambda = 0.3
    grad = x -> d .* x - b
    prox = (z, c) -> sign.(z) .* max.(abs.(z) .- lambda*c, 0)
    cost = x -> sum(d .* abs2.(x))/2 - sum(b .* x) + lambda*sum(abs, x)
    solution = prox(b, 1.) ./ d
    snapshot = (k, x, y, r) -> (k, copy(x), copy(y), r)
    for interval in (1, 7), rtol in (0., 1e-8)
        x, history, info = pogm_restart_gr(x0, cost, grad, 100.;
            g_prox=prox, niter=2000, atol=1e-7, rtol,
            check_every=interval, return_info=true, fun=snapshot)
        @test info.converged
        @test info.reason == :converged
        @test 0 < info.iterations < 2000
        @test info.iterations % interval == 0
        @test length(history) == info.iterations + 1
        @test x == history[end][2]
        @test x ≈ solution atol=2e-6
        residual = norm((x - prox(x - grad(x)/100, 1/100))*100)
        @test info.residual ≈ residual
        @test residual <= info.tolerance
        initial_residual = norm((x0 - prox(x0 - grad(x0)/100, 1/100))*100)
        @test info.tolerance ≈ 1e-7 + rtol*initial_residual
        # Compare with a longer-budget run: early exit does not recompute a last step.
        _, reference = pogm_restart_gr(x0, cost, grad, 100.;
            g_prox=prox, niter=2001, fun=snapshot)
        @test history == reference[1:length(history)]
    end

    # Relative-only stopping and convergence exactly at the final budget.
    _, _, relative = pogm_restart_gr(x0, cost, grad, 100.;
        g_prox=prox, niter=2000, rtol=1e-8, return_info=true)
    @test relative.converged && relative.residual <= relative.tolerance
    _, _, final = pogm_restart_gr([1.], x -> sum(abs2, x)/2, x -> x, 1.;
        niter=1, atol=0.6, return_info=true)
    @test final.converged && final.iterations == 1 && final.reason == :converged

    # At a composite optimum the smooth gradient need not be zero.
    @test norm(grad(solution)) > 0.1
    x, h, info = pogm_restart_gr(solution, cost, grad, 100.;
        g_prox=prox, atol=1e-12, return_info=true)
    @test info.converged && info.iterations == 0
    @test x == solution && length(h) == 1

    # Budget exhaustion still checks the special final step, even off cadence.
    x, h, info = pogm_restart_gr(x0, cost, grad, 100.;
        g_prox=prox, niter=2, atol=0., check_every=7, return_info=true)
    xr, hr = pogm_restart_gr(x0, cost, grad, 100.; g_prox=prox, niter=2)
    @test x == xr && h == hr
    @test !info.converged && info.reason == :niter && info.iterations == 2
    @test info.residual ≈ norm((x - prox(x - grad(x)/100, 1/100))*100)
    x, h, info = pogm_restart_gr(x0, cost, grad, 100.;
        g_prox=prox, niter=0, atol=0., return_info=true)
    @test x == x0 && length(h) == 1 && info.iterations == 0
    @test !info.converged && info.reason == :niter
    _, _, info = pogm_restart_gr(x0, cost, grad, 100.; niter=2, return_info=true)
    @test info.residual === nothing && info.tolerance === nothing

    for badgrad in (x -> fill(NaN, size(x)), x -> fill(Inf, size(x)))
        _, h, info = pogm_restart_gr(x0, cost, badgrad, 100.; atol=1e-6, return_info=true)
        @test !info.converged && info.reason == :nonfinite_residual
        @test info.iterations == 0 && length(h) == 1
    end
    # A nonfinite mapping after an update must not be called convergence either.
    calls = Ref(0)
    badlater = x -> (calls[] += 1; calls[] > 2 ? fill(NaN, size(x)) : grad(x))
    _, h, info = pogm_restart_gr(x0, cost, badlater, 100.; atol=0., return_info=true)
    @test !info.converged && info.reason == :nonfinite_residual
    @test info.iterations == 1 && length(h) == 2

    for options in ((; atol=-1.), (; atol=NaN), (; atol=Inf),
        (; rtol=-1.), (; rtol=Inf), (; check_every=0), (; niter=-1))
        @test_throws ArgumentError pogm_restart_gr(x0, cost, grad, 100.; options...)
    end
    @test_throws ArgumentError pogm_restart_gr(x0, cost, grad, 0.; atol=1e-6)

    # Small L1 logistic and squared-hinge problems with analytic optima.
    for (smooth, gradient, L, optimum) in (
        (x -> log1p(exp(-x[1])), x -> [-1/(1+exp(x[1]))], 0.25, log(4.)),
        (x -> max(0., 1-x[1])^2/2, x -> [-max(0., 1-x[1])], 1., 0.8),
    )
        soft = (z, c) -> sign.(z) .* max.(abs.(z) .- 0.2*c, 0)
        objective = x -> smooth(x) + 0.2*sum(abs, x)
        x, _, info = pogm_restart_gr([0.], objective, gradient, L;
            g_prox=soft, niter=1000, atol=1e-9, return_info=true)
        @test info.converged && info.iterations < 1000
        @test x[1] ≈ optimum atol=1e-7
        @test info.residual <= 1e-9
    end
end

include(joinpath(@__DIR__, "pogm_test_problems.jl"))

@testset "L1 classification reference accuracy" begin
    for kind in (:logistic, :squared_hinge)
        p = pogm_test_problem(kind)
        # Check the fixture's gradient independently of the solver.
        xcheck = collect(range(-0.2, 0.3; length=p.n))
        direction = collect(range(0.3, -0.1; length=p.n))
        step = 1e-6
        difference = (p.smooth(xcheck + step*direction) -
            p.smooth(xcheck - step*direction)) / (2*step)
        @test difference ≈ sum(p.grad(xcheck) .* direction) rtol=1e-6 atol=1e-9

        # Long fixed-budget original solver reference, not the new stopping rule.
        reference, _ = pogm_restart(p.x0, p.cost, p.grad, p.L;
            g_prox=p.prox, niter=5000, mom=:pogm, restart=:gr)
        @test p.residual(reference) < 1e-10
        reference_cost = p.cost(reference)
        for interval in (1, 5, 20)
            x, history, info = pogm_restart_gr(p.x0, p.cost, p.grad, p.L;
                g_prox=p.prox, niter=2000, atol=1e-7,
                check_every=interval, return_info=true)
            @test info.converged && info.reason == :converged
            @test info.iterations < 2000
            @test length(history) == info.iterations + 1
            @test p.residual(x) <= 1e-7
            @test info.residual ≈ p.residual(x)
            @test abs(p.cost(x) - reference_cost) < 1e-8
            @test norm(x - reference) / max(1., norm(reference)) < 1e-4
        end
    end
end

@testset "Fused updates preserve types and callback references" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64), shape in ((4,), (2, 2))
        x0 = reshape(T[1, -2, 3, -4], shape)
        original_x0 = copy(x0)
        gradient = x -> 0.25*x
        cost = x -> sum(abs2, x)/8
        prox = (z, c) -> z / (1 + 0.1*c)
        # Deliberately retain raw references: later iterations must not overwrite them.
        retain = (k, x, y, restarted) -> (k, x, y, restarted)
        x, history = pogm_restart_gr(x0, cost, gradient, 1.;
            g_prox=prox, niter=40, fun=retain)
        xr, reference = pogm_restart(x0, cost, gradient, 1.;
            g_prox=prox, niter=40, fun=retain)
        @test x == xr
        @test history == reference
        @test x0 == original_x0
        @test history[2][2] !== history[3][2]
        @test history[2][3] !== history[3][3]
    end
end

@testset "Optional mutating callbacks" begin
    snapshot = (k, x, y, r) -> (k, copy(x), copy(y), r)
    retain = (k, x, y, r) -> (k, x, y, r)
    unused = (args...) -> error("nonmutating callback should not be called")
    for kind in (:logistic, :squared_hinge)
        p = pogm_test_problem(kind)
        dest = similar(p.x0)
        p.grad!(dest, p.x0)
        @test dest ≈ p.grad(p.x0) rtol=1e-14
        p.prox!(dest, p.x0, 0.3)
        @test dest == p.prox(p.x0, 0.3)
        for stopping in (false, true)
            options = stopping ? (; atol=1e-7, check_every=5) : (;)
            xr, hr, ir = pogm_restart_gr(p.x0, p.cost, p.grad, p.L;
                g_prox=p.prox, niter=200, fun=snapshot, return_info=true, options...)
            for (use_grad, use_prox) in ((true, false), (false, true), (true, true))
                grad_calls, prox_calls = Ref(0), Ref(0)
                grad! = (d, x) -> (grad_calls[] += 1; p.grad!(d, x))
                prox! = (d, z, c) -> (prox_calls[] += 1; p.prox!(d, z, c))
                x, h, info = pogm_restart_gr(p.x0, p.cost, use_grad ? unused : p.grad, p.L;
                    g_prox=use_prox ? unused : p.prox,
                    f_grad! = use_grad ? grad! : nothing,
                    g_prox! = use_prox ? prox! : nothing,
                    niter=200, fun=retain, return_info=true, options...)
                @test x ≈ xr rtol=1e-12 atol=1e-14
                @test length(h) == length(hr)
                @test all(a[1] == b[1] && a[4] == b[4] &&
                    isapprox(a[2], b[2]; rtol=1e-12, atol=1e-14) &&
                    isapprox(a[3], b[3]; rtol=1e-12, atol=1e-14) for (a,b) in zip(h,hr))
                @test (info.iterations, info.reason) == (ir.iterations, ir.reason)
                @test grad_calls[] == (use_grad ? ir.iterations + (stopping ? 1+ir.iterations÷5 : 0) : 0)
                @test prox_calls[] == (use_prox ? ir.iterations + (stopping ? 1+ir.iterations÷5 : 0) : 0)
                @test h[2][2] !== h[3][2]
                @test h[2][3] !== h[3][3]
                @test p.x0 == zeros(p.n)
            end
        end
    end
    for T in (Float32, Float64, ComplexF32, ComplexF64), shape in ((4,), (2, 2))
        x0 = reshape(T[1, -2, 3, -4], shape)
        grad = x -> 0.25*x
        prox = (z, c) -> z / (1+0.1*c)
        grad! = (d, x) -> (d .= 0.25 .* x; nothing)
        prox! = (d, z, c) -> (d .= z ./ (1+0.1*c); nothing)
        cost = x -> sum(abs2, x)/8
        xr, hr = pogm_restart_gr(x0, cost, grad, 1.; g_prox=prox, niter=20, fun=snapshot)
        x, h = pogm_restart_gr(x0, cost, unused, 1.;
            f_grad! = grad!, g_prox! = prox!, niter=20, fun=retain)
        @test x == xr
        @test h == hr
    end
    # An optimum can be detected before the first iteration using only ! callbacks.
    x, h, info = pogm_restart_gr(zeros(3), x -> sum(abs2, x), unused, 1.;
        f_grad! = (d,x) -> fill!(d, 0), g_prox! = (d,z,c) -> copyto!(d,z),
        atol=0., return_info=true)
    @test info.converged && info.iterations == 0 && length(h) == 1
    @test x == zeros(3)
end
