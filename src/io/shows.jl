#=
shows.jl
2019-07-04, Jeff Fessler, University of Michigan
=#

export @shows

using Base: show_unquoted


function shows_string(file::String, line::String, varname::String, var::Any)
    size_str = typeof(var) <: Tuple ? "($(length(var)),)" : size(var)
    return "$file: $line $varname: $(typeof(var)) $(size_str)\n"
end

# print only if the session is interactive
function shows_print(
    string::String ;
    io::IO = isinteractive() ? stdout : IOBuffer(),
)
    print(io, string)
end


"""
    @shows expr

Show the type and size of an expression `expr` (typically a variable).

More concise output than `@show`
and typically this is all that is needed when debugging code.
"""
macro shows(ex)
    local line = __source__.line
    local file = basename(string(__source__.file))
    local blk = Expr(:block)
    push!(blk.args, :(shows_print(shows_string(
        string($(esc(file))),
        $(sprint(show_unquoted,line)),
        $(sprint(show_unquoted,ex)),
        $(esc(ex))
        ))))
#=
    push!(blk.args, :(print(string($(esc(file))) * ":")))
    push!(blk.args, :(print($(sprint(show_unquoted,line) * " "))))
    push!(blk.args, :(print($(sprint(show_unquoted,ex) * ": "))))
    push!(blk.args, :(print(string(typeof($(esc(ex)))))))
    push!(blk.args, :(print(" " * string(size($(esc(ex)))) * "\n")))
=#
    return blk
end


_unit_precision(::T) where {T <: Number} = T

"""
    _show_struct(io::IO, ::MIME"text/plain", st::Any)
Informative way to show fields of a struct (composite type).
"""
function _show_struct(io::IO, ::MIME"text/plain", st::Any)
    println(io, "$(typeof(st)) :")
    for f in fieldnames(typeof(st))
        p = getfield(st, f)
        t = (p isa Number) ? _unit_precision(p) : typeof(p)
        println(io, " ", f, "::", t, " ",
            p isa Number ? p :
            p isa AbstractArray ? size(p) :
            "",
        )
    end
end
