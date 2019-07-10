#=
plot/setup
2019-06-24, Jeff Fessler, University of Michigan
=#

# module MIRT_plot_setup # an experiment; would need Reexport

export mirt_setup_plot

using Plots: plot

"""
`function mirt_setup_plot()`

MIRT tests use Plots (with default GR backend).
When running on a remote linux server over X11 I get error messages from GKS.
The following kludge (?) tries to overcome this problem.
"""
function mirt_setup_plot()

	try
		plot(ones(3))
	catch
		if Sys.islinux
			ENV["GKS_WSTYPE"] = "x11"
		else
			@warn "Plot not working for Sys.KERNEL=$(Sys.KERNEL); test may fail"
		end
	end

	try
		plot(ones(3))
	catch
		@warn "Plot still not working; test may fail"
	end

	nothing
end

# end
