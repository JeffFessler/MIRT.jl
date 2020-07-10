# rmsd100
# root mean-squared difference between two arrays in %

export rmsd100

using LinearAlgebra: norm


"""
    rmsd = rmsd100(x, y ; mask)

Compute 100 * RMSD (root mean squared difference) between `x` and `y` within domain mask.

in
- `x` : array
- `y` : another array of same size

option:
- `mask::Array{Bool}` : domain over which to compute the RMSE; default `trues(size(x))`

out
- rmsd : rmsd of `x` vs `y` within `mask` in %
"""
function rmsd100(x::AbstractArray{<:Number}, y::AbstractArray{<:Number}
        ; mask::AbstractArray{Bool} = trues(size(x)))
    return 100 * norm(x[mask] - y[mask]) / sqrt(sum(mask))
end
