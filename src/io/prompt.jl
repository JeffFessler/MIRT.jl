#=
prompt.jl
prompt user to hit a key to continue
=#

export prompt

import Plots # gui, isplotnull

#using REPL #: RadioMenu, request
import REPL
using REPL.TerminalMenus

# global prompt state
prompt_state = :prompt


"""
`prompt()`
* prompt user to hit any key to continue, after `gui()`
* some keys have special actions: `[q]uit [d]raw [n]odraw`
* call `prompt(:prompt)` to revert to default
"""
function prompt( ; gui::Bool=true)
	global prompt_state

	prompt_state === :nodraw && return nothing

	gui && !Plots.isplotnull() && display(plot!()) # Plots.gui()

	(prompt_state === :draw) && return nothing
	!isinteractive() && return nothing

@show prompt_state

	c = wait_for_key()

	(c == 'd') && (prompt_state = :draw)
	(c == 'n') && (prompt_state = :nodraw)
	(c == 'q') && throw("quit")

@show prompt_state

	return nothing

#=
	list = ["continue*", "run", "draw", "quit", "nodraw"]
	menu = RadioMenu(list) #;, pagesize=3)

	preface = ""
	[name, line] = caller_name
	if isempty(name)
		preface = []
	else
		preface = sprintf('%s %d: ', name, line)
	end

	what = "enter to continue (or select)"

	choice = request(preface * what, menu)
	nothing
=#

end


"""
`function wait_for_key( ; prompt=?, io=stdin)`
from:
https://discourse.julialang.org/t/wait-for-a-keypress/20218
"""
function wait_for_key( ; io_in = stdin, io_out = stdout,
		prompt::String = "press any key [d]raw [n]odraw [q]uit : ")

	print(io_out, prompt)

	Base.Sys.iswindows() && (readline(); return nothing) # PC

	# non-windows version:
	setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), io_in.handle, raw)
	setraw!(true)
	char = Char(read(io_in, 1)[1])
	setraw!(false)
	write(io_out, char)
	write(io_out, "\n")

	return char
end


"""
`prompt(key::symbol)`
set prompt state to one of:
- `:prompt` call gui() if possible then prompt user
- `:draw` call gui() if possible then continue
- `:nodraw` do not call gui(), just continue

Actually it calls `display(plot!())` instead of `gui()`

`prompt(:test)`
self test
"""
function prompt(key::Symbol)
	global prompt_state

	if key === :state
		return prompt_state
	end

	if key === :test
		tmp = prompt(:state) # save current state
		prompt(:draw)
		@test prompt(:state) === :draw
		prompt()
		prompt(tmp) # return to original state
		return true
	end

	key ∉ (:prompt, :draw, :nodraw) && throw("prompt $prompt")
	prompt_state = key
end
