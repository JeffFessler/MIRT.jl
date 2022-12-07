# test/io/caller_name.jl

using MIRT: caller_name

using Test: @test


function caller_name_test()

    function f2()
        caller_name()
    end

    line = 2 + @__LINE__
    function f1()
        f2() # this is two lines below @__LINE__ above
    end

    @test isa(f1(), String)
    @test f1()[end-12:end] == ".jl $line f1(): "

    true
end

@test caller_name_test()
