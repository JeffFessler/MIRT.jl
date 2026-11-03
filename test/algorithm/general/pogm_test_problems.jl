# Deterministic synthetic fixtures shared by correctness tests and the standalone
# benchmark. No ECG data or external packages are needed.
using LinearAlgebra: opnorm, norm, mul!
using Random: MersenneTwister, randn

function pogm_test_problem(kind; m=240, n=40, seed=17)
    rng = MersenneTwister(seed)
    A = randn(rng, m, n) .* reshape(10. .^ range(0, -1; length=n), 1, n)
    truth = zeros(n)
    truth[1:5] = [2., -1.5, 1., -2., 0.5]
    labels = ifelse.(A * truth + 0.7 * randn(rng, m) .>= 0, 1., -1.)
    B = A .* labels
    lambda = 0.015
    # Stable log(1 + exp(z)) and sigmoid(-z).
    softplus = z -> max(z, 0.) + log1p(exp(-abs(z)))
    sigmoid_negative = z -> z >= 0 ? exp(-z)/(1+exp(-z)) : 1/(1+exp(z))
    if kind == :logistic
        smooth = x -> sum(softplus, -(B*x)) / m
        grad = x -> -(B' * sigmoid_negative.(B*x)) / m
        L = opnorm(B)^2 / (4*m)
    elseif kind == :squared_hinge
        smooth = x -> sum(abs2, max.(0., 1 .- B*x)) / (2*m)
        grad = x -> -(B' * max.(0., 1 .- B*x)) / m
        L = opnorm(B)^2 / m
    else
        throw(ArgumentError("unknown problem $kind"))
    end
    # Each fixture owns its scratch space; these callbacks are used sequentially.
    scores = zeros(m)
    grad! = function (dest, x)
        mul!(scores, B, x)
        if kind == :logistic
            scores .= sigmoid_negative.(scores)
        else
            @. scores = max(0., 1 - scores)
        end
        mul!(dest, B', scores)
        @. dest = -dest / m
        return nothing
    end
    prox! = function (dest, z, c)
        @. dest = sign(z) * max(abs(z) - lambda*c, 0.)
        return nothing
    end
    prox = (z, c) -> sign.(z) .* max.(abs.(z) .- lambda*c, 0.)
    cost = x -> smooth(x) + lambda*sum(abs, x)
    residual = x -> norm((x - prox(x - grad(x)/L, 1/L))*L)
    return (; x0=zeros(n), cost, grad, prox, L, residual, smooth, lambda, m, n, grad!, prox!)
end
