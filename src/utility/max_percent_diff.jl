#=
max_percent_diff.jl
Copyright 2000-9-16, Jeff Fessler, University of Michigan
=#

export max_percent_diff


"""
    max_percent_diff(s1, s2, [options])

Compute the "maximum percent difference" between two signals: `s1, s2`.

Default is to normalize by `maximum(abs, s1)`.

options
- `maxboth::Bool = false` use max of both arguments to normalize
- `normalize::Bool = false` normalize each before comparing
"""
function max_percent_diff(s1, s2 ;
    maxboth::Bool = false,
    normalize::Bool = false,
)

    # check to make sure that we have comparable signals
    axes(s1) == axes(s2) ||
        throw(DomainError("Dimension mismatch $(axes(s1)) $(axes(s2))"))

    any(isnan, s1) && @warn("NaN values in s1 input to max_percent_diff")
    any(isnan, s2) && @warn("NaN values in s2 input to max_percent_diff")

    if normalize
        normFactor = sum(abs, s2) / sum(abs, s1)
        isinteractive() &&
            (@info "Normalization factor (array2 / array1): $normFactor \n")
        iszero(sum(abs, s1)) && error("Need nonzero s1 norm")
        iszero(sum(abs, s2)) && error("Need nonzero s2 norm")
        s1 /= sum(abs, s1)
        s2 /= sum(abs, s2)
    end

    if maxboth
        denom = max(maximum(abs, s1), maximum(abs, s2))
        d = iszero(denom) ? 0 : maximum(abs, s1 - s2) / denom
    else
        denom = maximum(abs, s1)
        err = maximum(abs, s1 - s2)
        if iszero(denom)
            if iszero(err)
                d = 0
            else
                error("Cannot normalize with 0 array")
            end
        else
            d = err / denom
        end
    end

    return 100 * d
end
