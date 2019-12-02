"""
    exp_xform(x, u, v,mode = :matrix)
in:
	x	[N L]	complex vector(s)
	u	[D N]	complex vectors
	v	[D M]	complex vectors
   mode Symbol
out:
	y	[M L]	complex vector
	y(m,l) = sum_n x(n,l) exp(-sum(u(:,n) .* v(:,m)))
   Iterates through subsets of the ML matrix designated by :mode
   (i.e. row, column, element, or just computing the entire matrix)
	This is the 'slow' 'exact' transform model for MRI.
	All 3 inputs must be the same type (single or double).
	Output type will be same as input.
"""
function exp_xform(x::Array{<:Complex{<:Number}},u::Array{<:Complex{<:Number}},v::Array{<:Complex{<:Number}};mode::Symbol = :matrix)
   if(typeof(x) != typeof(u))
      throw("Types of first and second arguments do not match.")
   end
   if(typeof(x) != typeof(v))
      throw("Types of first and third arguments do not match.")
   end
   if(typeof(u) != typeof(v))
      throw("Types of second and third arguments do not match.")
   end
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
   else
      throw("Invalid mode parameter.")
   end
end
function exp_xform(x::Symbol)
   modes = [:element,:row,:column]
    if(x != :test)
       throw("Invalid argument for exp_xform.")
    end
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

    y1 = exp_xform(X,U,V,mode = :matrix)
    for i = 1:size(modes,1)
      print("Case $(modes[i])")
      @time y2 = exp_xform(X,U,V,mode = modes[i])
      d = max_percent_diff(y1, y2)
      print("double max % diff = $d\n")
      if d < 1e-12
        print("double appears to be working\n")
      else
        print("double may have a problem?\n")
      end

    if true
    	xs = Array{Complex{Float32},2}(X)
    	us = Array{Complex{Float32},2}(U)
    	vs = Array{Complex{Float32},2}(V)
        print("Single tests")
        @time y3 = exp_xform(X,U,V,mode = modes[i])
        d = max_percent_diff(y1, y3)
        print("single max % diff = $d\n")
        if d < 1e-4
          print("single appears to be working\n")
        else
          print("single may have a problem?\n")
        end
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
     print("Case $(modes[i])")
     @time double[:,:,i] = exp_xform(X,U,V,mode = modes[i])

   if true
     xs = Array{Complex{Float32},2}(X)
     us = Array{Complex{Float32},2}(U)
     vs = Array{Complex{Float32},2}(V)
      print("Single: ")
      @time single[:,:,i] = exp_xform(X,U,V,mode = modes[i])
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
  return true
end
