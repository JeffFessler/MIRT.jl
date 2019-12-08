#=
exp_xform.jl

2019-12-08 Connor Martin
Based on exp_xform_mex.c
Copyright 2004-9-23, Jeff Fessler, University of Michigan
=#

export exp_xform

using MIRT: max_percent_diff # todo
using Test: @test, @inferred

"""
    exp_xform(x, u, v ; mode::Symbol = :matrix)

in:
* `x [N L]` possibly complex vector(s)
* `u [D N]` possibly complex vectors
* `v [D M]` possibly complex vectors
* `mode::Symbol` `:matrix` (default) | `:element` | `:row` | `:column`

out:
* `y [M L]` typically complex vector
 `y[m,l] = sum_n x[n,l] exp(-sum(u[:,n] .* v[:,m]))`

Iterates through subsets of the ML matrix designated by `:mode`
(i.e. row, column, element, or just computing the entire matrix)
This is the 'slow' 'exact' transform model for MRI.

Output type will depend on input types.
"""
function exp_xform(x::AbstractMatrix{<:Number},
        u::AbstractMatrix{<:Number},
        v::AbstractMatrix{<:Number}
        ; mode::Symbol = :matrix)

    mode ∉ (:matrix, :element, :row, :column) && throw("Invalid mode parameter.")

    T = promote_type(eltype(u), eltype(v), eltype(x), ComplexF32)
@show T

    N = size(u,2)
    M = size(v,2)
    out = zeros(T, M, size(x,2))

    if (mode === :matrix)
        return exp.(-(transpose(v) * u)) * x # [M D] * [D N] * [N L]

    elseif (mode === :element)
        # avoid generating intermediate steps by iterating thru N and M
        # dot the columns
        for n = 1:N
            for m = 1:M
                t = transpose(u[:,n]) * v[:,m] #calculate one spot in the array
                out[m,:] .+= (exp.(-t) .* x[n,:])
            end
        end
        return out

    elseif (mode === :row) # iterate through rows of the large matrix
        for m = 1:M
            t = transpose(v[:,m]) * u # [1 D] * [D N] -> [1 N]
            # m == 1 && @show size(t), size(exp.(-t) * x)
            # product is a 2d matrix (row), but should be a column vector:
            out[m,:] .+= (exp.(-t) * x)[1,:]
        end
        return out

    elseif (mode === :column)
        for n = 1:N
            t = transpose(v) * u[:,n] # [M D] * [D .] -> [M]
            # n == 1 && @show size(t), size(x[n,:])
            out .+= (exp.(-t) * transpose(x[n,:])) # outer product [M 1] * [1 L]
        end
        return out

    end
end


# handle "vector" case where L=1
exp_xform(x::AbstractVector{<:Number},
        u::AbstractMatrix{<:Number},
        v::AbstractMatrix{<:Number}
        ; kwargs...) = exp_xform(reshape(x, :, 1), u, v ; kwargs...)[:,1]


"""
    exp_xform(:test)
self test
"""
function exp_xform(x::Symbol ; time::Bool = false)
    x != :test && throw("Invalid argument for exp_xform.")

    modes = [:element, :row, :column]
    N = 500
    M = 6000
    D = 3
    L = 1
    X = Complex.(randn(N,L), randn(N,L))
    U = Complex.(randn(D,N), randn(D,N))
    V = Complex.(randn(D,M), randn(D,M))
    # todo: cut this and test F32 and F64 separately in loop
    if true
    	X = Array{Complex{Float64},2}(Array{Complex{Float32},2}(X))
    	U = Array{Complex{Float64},2}(Array{Complex{Float32},2}(U))
    	V = Array{Complex{Float64},2}(Array{Complex{Float32},2}(V))
    end

    # todo: @inferred
    y1 = exp_xform(X,U,V,mode = :matrix)
    for i = 1:size(modes,1)
      print("Case $(modes[i])")
      time && @time y2 = exp_xform(X,U,V,mode = modes[i])
      !time && (y2 = exp_xform(X,U,V,mode = modes[i]))
      d = max_percent_diff(y1, y2)
      print("double max % diff = $d\n")
      d < 1e-12 && print("double appears to be working\n")
      d >= 1e-12 && print("double may have a problem?\n")
    if true
    	xs = Array{Complex{Float32},2}(X)
    	us = Array{Complex{Float32},2}(U)
    	vs = Array{Complex{Float32},2}(V)
        print("Single tests")
        time && @time y3 = exp_xform(X,U,V,mode = modes[i])
        !time && (y3 = exp_xform(X,U,V,mode = modes[i]))
        d = max_percent_diff(y1, y3)
        print("single max % diff = $d\n")
        d >= 1e-4 && print("single may have a problem?\n")
        d < 1e-4 && print("single appears to be working\n")
      end
      print("\n\n")
   end

   print("stress tests (not really) with medium N/M coeff\n")
   N = 10000
   M = 100
   D = 20
   L = 10
   X = Complex.(randn(N,L),randn(N,L))
   U = Complex.(randn(D,N),randn(D,N))
   V = Complex.(randn(D,M),randn(D,M))
   single = Array{Complex{Float32}}(zeros(M,L,size(modes,1)))
   double = Array{Complex{Float64}}(zeros(M,L,size(modes,1)))
   for i = 1:size(modes,1)
     time && print("Case $(modes[i]) \n")
     time && @time double[:,:,i] = exp_xform(X,U,V,mode = modes[i])
     !time && (double[:,:,i] = exp_xform(X,U,V,mode = modes[i]))
     if true
     xs = Array{Complex{Float32},2}(X)
     us = Array{Complex{Float32},2}(U)
     vs = Array{Complex{Float32},2}(V)
      time && print("Single: ")
      time && @time single[:,:,i] = exp_xform(X,U,V,mode = modes[i])
      !time && (single[:,:,i] = exp_xform(X,U,V,mode = modes[i]))
     end
     print("\n\n")
  end
  for i = 1:size(modes,1)
     for j = 1:(i-1)
        d = max_percent_diff(double[:,:,i],double[:,:,j])
        print("double max % diff between $(modes[i]) and $(modes[j]) = $d\n")
        d = max_percent_diff(single[:,:,i],single[:,:,j])
        print("single max % diff between $(modes[i]) and $(modes[j]) = $d\n")
     end
  end
  true
end

exp_xform(:test)
