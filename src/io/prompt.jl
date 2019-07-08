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
prompt user to hit enter to continue, after gui()
"""
function prompt(; gui::Bool=true)
	global prompt_state

	prompt_state == :nodraw && return nothing

	gui && !Plots.isplotnull() && Plots.gui()

	prompt_state == :prompt && wait_for_key()
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


"""
`prompt(key::symbol)`
set prompt state to one of:
- `:prompt` call gui() if possible then prompt user
- `:draw` call gui() if possible then continue
- `:nodraw` do not call gui(), just continue

`prompt(:test)`
self test
"""
function prompt(key::Symbol)
	global prompt_state

	if key == :state
		return prompt_state
	end

	if key == :test
		prompt(:draw)
		@test prompt(:state) == :draw
		prompt()
		return true
	end

	key ∉ (:prompt, :draw, :nodraw) && throw("prompt $prompt")
	prompt_state = key
end
