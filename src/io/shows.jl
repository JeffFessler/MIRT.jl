#=
shows.jl
2019-07-04, Jeff Fessler, University of Michigan
=#

export @shows

using Base: show_unquoted


function shows_string(file::String, line::String, varname::String, var::Any)
	return "$file: $line $varname: $(typeof(var)) $(size(var))\n"
end

# print only if the session is interactive
function shows_print(string::String
		; io::IO = isinteractive() ? stdout : IOBuffer())
	print(io, string)
end


"""
`@shows expr`

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


#=
This test was used solely during development

"""
`shows(:test)
test macro `@shows`
"""
function shows(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))
	var = ones(3,4)
	@shows var
	true
end
=#

# shows(:test)
