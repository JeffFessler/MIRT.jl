# prompt.jl
# prompt user to hit a key to continue

import Plots # gui

#using REPL #: RadioMenu, request
import REPL
using REPL.TerminalMenus

"""
`prompt()`
prompt user to hit enter to continue, after gui()
"""
function prompt(; gui::Bool=true)
	gui && Plots.gui()
	wait_for_key()
	return nothing

	list = ["continue*", "run", "draw", "quit", "nodraw"]
	menu = RadioMenu(list) #;, pagesize=3)

	preface = ""
#=
	[name, line] = caller_name
	if isempty(name)
		preface = []
	else
		preface = sprintf('%s %d: ', name, line)
	end
=#

	what = "enter to continue (or select)"

	choice = request(preface * what, menu)
	nothing
end


"""
`function wait_for_key(; prompt=?, io=stdin)`
from:
https://discourse.julialang.org/t/wait-for-a-keypress/20218
"""
function wait_for_key(; prompt::String = "press any key: ", io = stdin)
	setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), io.handle, raw)
	print(io, prompt)
	setraw!(true)
	read(io, 1)
	setraw!(false)
	nothing
end
