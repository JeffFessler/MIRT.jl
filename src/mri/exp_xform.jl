#=
exp_xform.jl
=#
export exp_xform

"""
    exp_xform(x, u, v, mode::Symbol = :matrix)
in:
	`x	[N L]`	complex vector(s)
	`u	[D N]`	complex vectors
	`v	[D M]`	complex vectors
   `mode::Symbol` `:matrix` (default) | `:element` | `:row` | `:column`
out:
	`y	[M L]`	complex vector
	`y[m,l] = sum_n x[n,l] exp(-sum(u[:,n] .* v[:,m]))`
   Iterates through subsets of the ML matrix designated by `:mode`
   (i.e. row, column, element, or just computing the entire matrix)
	This is the 'slow' 'exact' transform model for MRI.
	All 3 inputs must be the same type (single or double).
	Output type will be same as input.
"""
function exp_xform(x::AbstractArray{<:Number},
        u::AbstractArray{<:Number},
        v::AbstractArray{<:Number}
        ; mode::Symbol = :matrix)
   mode != :matrix && mode != :element && mode != :row && mode != :column && throw("Invalid mode parameter.")
   out = Array{promote_type(ComplexF32,typeof(x))}(zeros(size(v,2),size(x,2)))

   if(mode === :matrix)
        return exp.(-1 .* (transpose(v) * u)) * x

   elseif(mode === :element)
      #avoid generating intermediate steps by iterating thru N and M
      #dot the columns
      #@show size(x)
      for i = 1:size(u,2)
        for j = 1:size(v,2)
           t = transpose(u[:,i]) * v[:,j] #calculate one spot in the array
           p = out[j,:] .+ (exp.(-1 * t) .* x[i,:])
           out[j,:] .= p
        end
     end
      return out

   elseif(mode === :row)
      #iterate through the rows of the large matrix
      #@show size(x)
     for j = 1:size(v,2)
           t = transpose(v[:,j]) * u
           p = out[j,:] .+ (exp.(-1 * t) * x)[1,:] #this output is a 2d matrix, but should be a column vector (it's a row)
           out[j,:] .= p
     end
     return out

  elseif(mode === :column)
      for i = 1:size(u,2)
            t = transpose(v) * u[:,i]
            #@show size(exp.(-1 * t))
            #@show size(x[i,:])
            p = out .+ (exp.(-1 * t) * transpose(x[i,:]))#outer product
            out = p
      end
      return out
   end
end


"""
    exp_xsform(:test)
self test
"""
function exp_xform(x::Symbol; time = false)
   modes = [:element,:row,:column]
    x != :test && throw("Invalid argument for exp_xform.")
    N = 500
    M = 6000
    D = 3
    L = 1
    X = Complex.(randn(N,L), randn(N,L))
    U = Complex.(randn(D,N), randn(D,N))
    V = Complex.(randn(D,M), randn(D,M))
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
