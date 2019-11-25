#=
max_percent_diff.jl
Copyright 2000-9-16, Jeff Fessler, University of Michigan
=#

export max_percent_diff

using Test: @test


"""
    d = max_percent_diff(s1, s2, [options])

compute the "maximum percent difference" between two signals: `s1, s2`
default is to normalize by `maximum(abs.(s1))`

options
- `maxboth::Bool` use max of both arguments to normalize; default `false`
- `normalize::Bool` normalize each before comparing; default `false`

"""
function max_percent_diff(s1, s2 ;
        maxboth::Bool = false, normalize::Bool = false)

    #check to make sure that we have comparable signals
    size(s1) != size(s2) &&
        throw(DomainError("Dimension mismatch $(size(s1)) $(size(s2))"))

    any(isnan.(s1)) && @warn("NaN values in s1 input to max_percent_diff")
    any(isnan.(s2)) && @warn("NaN values in s2 input to max_percent_diff")

    if normalize
        normFactor = sum(abs.(s2)) / sum(abs.(s1))
        print("Normalization factor (matrix 2 / matrix 1): $normFactor \n")
        s1 /= sum(abs.(s1))
        s2 /= sum(abs.(s2))
    end

    if maxboth
        denom = max(maximum(abs.(s1)), maximum(abs.(s2)))
        d = denom == 0 ? 0 : maximum(abs.(s1 - s2)) / denom
    else
        denom = maximum(abs.(s1))
        if denom == 0
            denom = maximum(abs.(s2))
        end
        d = denom == 0 ? 0 : maximum(abs.(s1 - s2)) / denom
    end

    return 100 * d
end


function max_percent_diff_test()
    @test max_percent_diff([200], [202]) == 1
    @test max_percent_diff([0 200], [0 202]) == 1
    @test max_percent_diff([0 200], [0 198], maxboth=true) == 1
    @test max_percent_diff([0 100], [0 200], normalize=true) == 0
    true
end


"""
    max_percent_diff(:test)
self test
"""
function max_percent_diff(test::Symbol)
    test != :test && throw(DomainError("Invalid signal"));
    return max_percent_diff_test()
end

# max_percent_diff(:test)
