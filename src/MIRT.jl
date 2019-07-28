# MIRT.jl


"""
`MIRT` is the Michigan Image Reconstruction Toolbox
"""
module MIRT
	using Reexport

	@reexport using MIRTio # make I/O routines available

	include("z-all.jl")

end # module
