# test/io/shows.jl

using MIRT: @shows, _show_struct
using Test: @test

var = ones(3,4)
@test (@shows var) isa Nothing

struct Tester
    a
    b
    c
end

Base.show(io::IO, ::MIME"text/plain", src::Tester) =
   _show_struct(io, MIME("text/plain"), src)
tmp = Tester(5, 1:3, :test)
show(isinteractive() ? stdout : devnull, MIME("text/plain"), tmp)
