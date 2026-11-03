# POGM + gradient restart, isolated from pogm_restart (Kim and Fessler).

export pogm_restart_gr

"""
    x, out = pogm_restart_gr(x0, Fcost, f_grad, f_L; kwargs...)

POGM with gradient restart only. This is the isolated `mom=:pogm,
restart=:gr` path of [`pogm_restart`](@ref), with the same options except
`mom` and `restart`. Returns the secondary iterate and callback history.

By default runs exactly `niter` iterations, including the special final momentum step
when `f_mu == 0`. The callback receives `(iter, xk, yk, is_restart)` with
secondary `xk` and primary `yk`; its return value is saved, not used to stop.
`Fcost` is still evaluated initially and at each iteration to preserve the
original call behavior, although gradient restart does not use its value.
For compatibility, `niter=0` returns the original empty vector when stopping
is disabled; with stopping enabled, it checks and returns `x0` instead.

Set `atol` (absolute tolerance) or positive `rtol` to enable early stopping.
With `alpha = 1/f_L`, checks use the standard proximal-gradient mapping
`G(x) = (x - g_prox(x - alpha*f_grad(x), alpha)) / alpha` at the
returned secondary iterate, requiring
`norm(G(x)) <= atol + rtol*norm(G(x0))` (`atol=nothing` means zero).
The initial point is checked, then every `check_every` iterations and at
`niter`. Checks require extra gradient/prox evaluations. Gradient/prox
functions must return their results without mutating their inputs.
Early stopping returns the checked iterate without recomputing a special
final step; that step remains tied to the iteration budget `niter`.

With `return_info=true`, returns `(x, out, info)`, where `info` contains
`iterations`, `converged`, `reason` (`:converged`, `:niter`, or
`:nonfinite_residual`), `residual`, and `tolerance`. Residual and tolerance
are `nothing` when checking is disabled. History is truncated on stopping.
Elementwise solver updates use fused broadcasts to reduce temporary arrays.
Iterate arrays are still newly allocated, so callback references are not reused
as working buffers. Optionally provide `f_grad!(dest, x)` and/or `g_prox!(dest, z, c)`.
Each must fill `dest` completely without modifying its inputs; return values
are ignored. A supplied mutating callback replaces its nonmutating counterpart.
The gradient buffer has the shape of `x0` and element type
`promote_type(eltype(x0), typeof(1. / f_L))`; its contents are overwritten on
subsequent gradient evaluations. These callbacks must support that working
precision and must not retain or alias solver buffers. The proximal destination
is fresh each time, distinct from `z`, preserving callback history references.
The existing positional `f_grad` argument remains required but is not called
when `f_grad!` is supplied. User callbacks may own additional scratch arrays.
This reduces callback allocations; the solver is still not allocation-free.
"""
function pogm_restart_gr(
    x0,
    Fcost::Function,
    f_grad::G,
    f_L::Real ;
    f_mu::Real = 0.,
    restart_cutoff::Real = 0.,
    bsig::Real = 1,
    niter::Int = 10,
    g_prox::P = (z, c::Real) -> z,
    f_grad!::Union{Nothing,Function} = nothing,
    g_prox!::Union{Nothing,Function} = nothing,
    atol::Union{Nothing,Real} = nothing,
    rtol::Real = 0.,
    check_every::Int = 1,
    return_info::Bool = false,
    fun::Function = (iter::Int, xk, yk, is_restart::Bool) -> undef,
) where {G<:Function,P<:Function}
    f_L < 0 && throw(ArgumentError("f_L=$f_L < 0"))
    f_mu < 0 && throw(ArgumentError("f_mu=$f_mu < 0"))
    bsig < 0 && throw(ArgumentError("bsig=$bsig < 0"))
    !((-1 < restart_cutoff) && (restart_cutoff < 1)) &&
        throw(ArgumentError("restart_cutoff=$restart_cutoff"))

    niter < 0 && throw(ArgumentError("niter=$niter < 0"))
    check_every < 1 && throw(ArgumentError("check_every must be positive"))
    atol !== nothing && !(isfinite(atol) && atol >= 0) &&
        throw(ArgumentError("atol must be finite and nonnegative"))
    !(isfinite(rtol) && rtol >= 0) &&
        throw(ArgumentError("rtol must be finite and nonnegative"))
    stopping = atol !== nothing || rtol > 0
    stopping && !(isfinite(f_L) && f_L > 0) &&
        throw(ArgumentError("early stopping requires finite f_L > 0"))

    # Optional mutating callbacks are adapted to the same algorithm below.
    # Only the gradient buffer is reused: proximal results can escape via fun.
    gradient = if f_grad! === nothing
        f_grad
    else
        T = promote_type(eltype(x0), typeof(1. / f_L))
        buffer = similar(x0, T)
        x -> (f_grad!(buffer, x); buffer)
    end
    proximal = if g_prox! === nothing
        g_prox
    else
        (z, c) -> (dest = similar(z); g_prox!(dest, z, c); dest)
    end

    L = f_L
    mu = f_mu
    q = mu/L
    told = 1
    sig = 1
    zetaold = 1 # dummy
    xold = x0
    yold = x0
    uold = x0
    zold = x0
    Fcost(x0)
    Fgradold = zeros(size(x0)) # dummy

    out = Array{Any}(undef, niter+1)
    out[1] = fun(0, x0, x0, false)
    xnew = stopping ? x0 : [] # preserve the original zero-budget behavior by default
    iterations = 0
    converged = false
    reason = :niter
    residual = nothing
    tolerance = nothing
    # This mapping is deliberately separate from the restart mapping below.
    mapping_norm = x -> norm((x .- proximal(x .- (1. / L) .* gradient(x), 1. / L)) .* L)
    if stopping
        residual = mapping_norm(x0)
        tolerance = (atol === nothing ? 0. : atol) +
            (iszero(rtol) ? 0. : rtol * residual)
        if !isfinite(residual)
            reason = :nonfinite_residual
        elseif residual <= tolerance
            converged = true
            reason = :converged
        end
        if reason != :niter
            resize!(out, 1)
            info = (; iterations, converged, reason, residual, tolerance)
            return return_info ? (xnew, out, info) : (xnew, out)
        end
    end

    for iter in 1:niter
        alpha = 1. / L
        fgrad = gradient(xold)
        is_restart = false
        unew = @. xold - alpha * fgrad

        if mu != 0
            beta = (2 + q - sqrt(q^2+8*q))^2 / 4. / (1-q)
        else
            if iter == niter
                tnew = 0.5 * (1 + sqrt(1 + 8 * told^2))
            else
                tnew = 0.5 * (1 + sqrt(1 + 4 * told^2))
            end
            beta = (told - 1) / tnew
        end

        if mu != 0
            gamma = (2 + q - sqrt(q^2+8*q)) / 2.
        else
            gamma = sig * told / tnew
        end

        znew = @. (unew + beta * (unew - uold) + gamma * (unew - xold)
                - beta * alpha / zetaold * (xold - zold))
        zetanew = alpha * (1 + beta + gamma)
        xnew = proximal(znew, zetanew)

        # Preserve the original non-standard POGM mapping (not G_alpha(xnew)).
        Fgrad = @. fgrad - 1/zetanew * (xnew - znew)
        ynew = @. xold - alpha * Fgrad
        Fcost(xnew)

        if gr_restart(Fgrad, ynew-yold, restart_cutoff)
            tnew = 1
            sig = 1
            is_restart = true
        elseif sum(Float64, real.(Fgrad .* Fgradold)) < 0
            sig = bsig * sig
        end
        Fgradold = Fgrad
        uold = unew
        zold = znew
        zetaold = zetanew

        out[iter+1] = fun(iter, xnew, ynew, is_restart)
        iterations = iter
        if stopping && (iter % check_every == 0 || iter == niter)
            residual = mapping_norm(xnew)
            if !isfinite(residual)
                reason = :nonfinite_residual
            elseif residual <= tolerance
                converged = true
                reason = :converged
            end
            if reason != :niter
                resize!(out, iter+1)
                break
            end
        end
        xold = xnew
        yold = ynew
        if iszero(mu)
            told = tnew
        end
    end

    info = (; iterations, converged, reason, residual, tolerance)
    return return_info ? (xnew, out, info) : (xnew, out)
end
