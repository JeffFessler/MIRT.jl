# test/mri/exp_xform.jl

using MIRT: exp_xform
using Test: @test, @test_throws, @inferred


# test for given data type
function exp_xform_test( ; T::Type{<:Complex{<:AbstractFloat}} = ComplexF32)

    modes = (:element, :row, :column)

    isinteractive() && (@info "1D tests with $T")
    N = 500
    M = 6000
    D = 3
    X = randn(T, N)
    U = randn(T, D, N)
    V = randn(T, D, M)
    y1 = @inferred exp_xform(X, U, V ; mode = :matrix)

    for mode in modes
        y2 = @inferred exp_xform(X, U, V ; mode=mode)
        @test y1 ≈ y2
    end

    isinteractive() && (@info "2D tests with $T")
    N = 10000
    M = 100
    D = 20
    L = 10
    X = randn(T, N, L)
    U = randn(T, D, N)
    V = randn(T, D, M)
    y1 = @inferred exp_xform(X, U, V ; mode = :matrix)

    for mode in modes
        y2 = @inferred exp_xform(X, U, V ; mode=mode)
        @test y1 ≈ y2
    end

    true
end


@test_throws String exp_xform(ones(2,2), ones(2,2), ones(2,2) ; mode=:bad)

# for T in (ComplexF32, ComplexF64)
    @test exp_xform_test() # T=T
# end
