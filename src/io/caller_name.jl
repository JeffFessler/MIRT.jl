#=
caller_name.jl
Use Stacktrace to show information about the
file, line number and calling function

Useful for messages and debugging.

Jeff Fessler, 2019-06-27, University of Michigan
=#

export caller_name


"""
`caller_name() or caller_name(;level=4)`

Return "filename line fun():" as a string
to describe where this function was called.

Stack levels:
- 1: #caller_name
- 2: caller_name()
- 3: function that invoked caller()
- 4: the calling function we want to return
Hence the default level is 4,
but we increment it by one in case user says `@show caller_name()`
in which case stack[3] is a macro expansion.
"""
function caller_name( ; level::Int = 4)
	stack = stacktrace()
	level += (stack[3].func == Symbol("macro expansion")) # trick
	(level < 1 || level > length(stack)) && throw("bad level $level")
	frame = stack[level]
	file = string(frame.file)
	file = basename(file)
	line = frame.line
	func = frame.func
#	display("$which, $file $line $func(): ")

	return "$file $line $func(): "
end
