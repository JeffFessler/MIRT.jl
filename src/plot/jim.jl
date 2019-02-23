# jim.jl
using Plots
using MosaicViews

"""
`jim(x, ...)`

jiffy image display of `x` using `heatmap`

in
* `x` image, can be 2D or higher, if higher then it uses `mosaicviews`

option
* `aspect_ratio` for heatmap; default `:equal`
* `clim` for heatmap; default `(minimum(x),maximum(x))`
* `color` colormap; default `:grays`
* `ncol` for mosaicview for 3D and higher arrays; default `0` does auto select
* `padval` padding value for mosaic view; default `(minimum(x)`
* `title` for heatmap; default `""`
* `xlabel` for heatmap; default `""`
* `ylabel` for heatmap; default `""`
* `yflip` for heatmap; default `true`
* `xtick` for heatmap; default `[1,size(x,1)]`
* `ytick` for heatmap; default `[1,size(x,2)]`

out
* returns plot handle

2019-02-23 Jeff Fessler, University of Michigan
"""
function jim(x;
	aspect_ratio = :equal,
	clim = [],
	color = :grays,
	ncol = 0,
	padval = minimum(x),
	title = "",
	xlabel = "",
	ylabel = "",
	xtick = [1, size(x,1)],
	ytick = [1, size(x,2)],
	yflip = true)

	if !isreal(x)
		@warn "magnitude"
		x = abs.(x)
	end

	if isempty(clim) # must wait until after possible abs() to do this
		clim = (minimum(x), maximum(x))
	end

	if ndims(x) > 2
		if ncol == 0
			ncol = Int(floor(sqrt(prod(size(x)[3:end]))))
		end
		z = mosaicview(x, padval, ncol=ncol, npad=1)
	else
		z = x
	end

heatmap(z', transpose=false,
	aspect_ratio=aspect_ratio,
	clim=clim,
	color=color,
	title=title,
	yflip=true,
	xlabel=xlabel,
	ylabel=ylabel,
	xtick=xtick,
	ytick=ytick)

end # jim


# show docstring if user calls it with no arguments
function jim()
	@doc jim
end
