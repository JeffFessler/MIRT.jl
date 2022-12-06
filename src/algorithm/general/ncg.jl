#=
ncg.jl
Nonlinear CG optimization
2019-03-16, Jeff Fessler, University of Michigan
=#

export ncg

using LinearAlgebra: I, norm, dot

# todo: ncg!


"""
`(x,out) = ncg(B, gradf, curvf, x0 ; ...)`

Nonlinear preconditioned conjugate gradient algorithm
to minimize a general "inverse problem" cost function of the form
``\\Psi(x) = \\sum_{j=1}^J f_j(B_j x)``
where each function ``f_j(v)`` has a quadratic majorizer of the form
```math
q_j(v;u) = f_j(u) + \\nabla f_j(u) (v - u) + 1/2 \\|v - u\\|^2_{C_j(u)}
```
where ``C_j(u)`` is diagonal matrix of curvatures.
(It suffices for each ``f_j`` to have a Lipschitz smooth gradient.)

This CG method uses a majorize-minimize (MM) line search.

in
- `B`		vector of ``J`` blocks ``B_1,...,B_J``
- `gradf`	vector of ``J`` functions return gradients of ``f_1,...,f_J``
- `curvf`	vector of ``J`` functions `z -> curv(z)` that return a scalar
or a vector of curvature values for each element of ``z``
- `x0`	initial guess; need `length(x) == size(B[j],2)` for ``j=1...J``

Usually `x0` is a `Vector` but it can be an `Array`
if each `B_j` is a linear operator (e.g., `LinearMapAO`)
of suitable "dimensions".

option
- `niter`	# number of outer iterations; default 50
- `ninner`	# number of inner iterations of MM line search; default 5
- `P`		# preconditioner; default `I`
- `betahow`	"beta" method for the search direction; default `:dai_yuan`
- `fun`		User-defined function to be evaluated with two arguments (x,iter).
   * It is evaluated at (x0,0) and then after each iteration.

output
- `x`		final iterate
- `out`		`[niter+1] (fun(x0,0), fun(x1,1), ..., fun(x_niter,niter))`
   * (all 0 by default). This is an array of length `niter+1`
"""
function ncg(
    B::AbstractVector{<:Any},
    gradf::AbstractVector{<:Function},
    curvf::AbstractVector{<:Function},
    x0::AbstractArray{<:Number} ; # usually a Vector
    niter::Int = 50,
    ninner::Int = 5,
    P = I, # trick: this is an overloaded I (by LinearMapsAA)
    betahow::Symbol = :dai_yuan,
    fun::Function = (x,iter) -> 0,
)

    Base.require_one_based_indexing(B, gradf, curvf)

	out = Array{Any}(undef, niter+1)
	out[1] = fun(x0, 0)

	J = length(B)

	x = copy(x0)
	dir = []
	grad_old = []
	grad_new = []

	Bx = [B[j] * x for j in 1:J] # u_j in course notes
	grad = (Bx) -> sum([B[j]' * gradf[j](Bx[j]) for j in 1:J])

for iter in 1:niter
	grad_new = grad(Bx) # gradient
	npgrad = -(P * grad_new)
	if iter == 1
		dir = npgrad
	else
		if betahow === :dai_yuan
			denom =	dot(grad_new - grad_old, dir)
			if denom == 0
				betaval = 0
			else
				betaval = dot(grad_new, P * grad_new) / denom
			end
		else
			throw(ArgumentError("unknown beta choice: $betahow"))
		end
		dir = npgrad + betaval * dir # search direction
	end
	grad_old = grad_new

	# MM-based line search for step size alpha
	# using h(a) = sum_j f_j(uj + a vj)
	Bd = [B[j] * dir for j in 1:J] # v_j in course notes

	alf = 0
	for ii in 1:ninner
	#	derh = alf -> sum([Bd[j]' * gradf[j](Bx[j] + alf * Bd[j]) for j in 1:J])
		derh = 0 # derivative of h(a)
		curv = 0
		for j in 1:J
			tmp = Bx[j] + alf * Bd[j]
			derh += real(dot(Bd[j], gradf[j](tmp)))
			curv += sum(curvf[j](tmp) .* abs2.(Bd[j]))
		end
		curv < 0 && throw("bug: curv < 0")
		if curv > 0
			alf = alf - derh / curv
		end
		if alf == 0
			break
		end
	end

	x += alf * dir
	for j in 1:J # update Bj * x
		Bx[j] += alf * Bd[j]
	end
	out[iter+1] = fun(x, iter)
end

	return x, out
#	return eltype(x0).(x), out # todo
end


"""
`(x,out) = ncg(grad, curv, x0, ...)`

special case of `ncg` (nonlinear CG) for minimizing a cost function
whose gradient is `grad(x)`
and that has a quadratic majorizer with diagonal Hessian given by
`curv(x)`.
Typically `curv = (x) -> L` where `L` is the Lipschitz constant of `grad`
"""
function ncg(
    grad::Function,
    curv::Function,
    x0::AbstractArray{<:Number} ;
    kwargs...,
)

	return ncg([I], [grad], [curv], x0; kwargs...)
end
