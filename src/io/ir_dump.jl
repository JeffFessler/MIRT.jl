#=
ir_dump.jl
2019-03-05 Jeff Fessler
=#

export ir_dump


"""
    ir_dump(x::Any ; io::IO = stdout)
    ir_dump(io::IO, x::Any)

Show all the fields of a structure or `NamedTuple` more nicely than dump() does
"""
function ir_dump(y::Any ; io::IO = stdout, ntuplemax::Int = 3)
    x = typeof(y)
    print(io, x)
    fields = fieldnames(x)
    fieldtypes = x.types
    for (idx,fd) in enumerate(fields)
        println(io)
        print(io, " ", fd, "::")
        ft = fieldtypes[idx]
        print(io, ft)
        if ft <: Number
            print(io, " ", getfield(y, fd))
        end
        if ft == String
            print(io, " \"", getfield(y, fd), "\"")
        end
        if ft <: Tuple # nice for ImageGeom
            tmp = getfield(y, fd)
            if length(tmp) <= ntuplemax
                print(io, " ", tmp)
            else
                print(io, " (")
                for it = 1:ntuplemax
                    print(io, tmp[it], ", ")
                end
                print(io, "…)")
            end
        end
        if ft <: AbstractArray{Bool} # nice for ImageGeom
            tmp = getfield(y, fd)
            print(io, " {", sum(tmp), " of ", length(tmp), "}")
        end
    end
    println(io)
    nothing
end


ir_dump(io::IO, x::Any) = ir_dump(x ; io=io)
