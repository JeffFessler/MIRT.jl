# jim.jl

using Plots
using MosaicViews

"""
`jim(z, ...)`

jiffy image display of `x` using `heatmap`

in
* `z` image, can be 2D or higher, if higher then it uses `mosaicviews`

option
* `aspect_ratio` for heatmap; default `:equal`
* `clim` for heatmap; default `(minimum(z),maximum(z))`
* `color` colormap; default `:grays`
* `ncol` for mosaicview for 3D and higher arrays; default `0` does auto select
* `padval` padding value for mosaic view; default `(minimum(z)`
* `title` for heatmap; default `""`
* `xlabel` for heatmap; default `""`
* `ylabel` for heatmap; default `""`
* `yflip` for heatmap; default `true`
* `x` for x axis; default `1:size(z,1)`
* `y` for y axis; default `1:size(z,2)`
* `xtick` for heatmap; default `[minimum(x),maximum(x)]`
* `ytick` for heatmap; default `[minimum(y),maximum(y)]`

out
* returns plot handle

2019-02-23 Jeff Fessler, University of Michigan
"""
function jim(z;
	aspect_ratio = :equal,
	clim = [],
	color = :grays,
	ncol = 0,
	padval = [],
	title = "",
	xlabel = "",
	ylabel = "",
	x = (1:size(z,1)),
	y = (1:size(z,1)),
	xtick = [minimum(x),maximum(x)],
	ytick = [minimum(y),maximum(y)],
	yflip = true)

	if !isreal(z)
		@warn "magnitude"
		z = abs.(z)
	end

	if isempty(clim) # must wait until after possible abs() to do this
		clim = (minimum(z), maximum(z))
	end

	if isempty(padval) # must wait until after possible abs() to do this
		padval = minimum(z)
	end

	if ndims(z) > 2
		if ncol == 0
			ncol = Int(floor(sqrt(prod(size(z)[3:end]))))
		end
		z = mosaicview(z, padval, ncol=ncol, npad=1)
	end

heatmap(x, y, z', transpose=false,
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



"""
`jim(x, y, z; kwargs...)`
"""
function jim(x, y, z; kwargs...)
	return jim(z; x=x, y=y, kwargs...)
end


# show docstring if user calls it with no arguments
function jim()
	@doc jim
end


function test_jim()
	jim(ones(4,3), title="test2")
	jim(ones(4,3,5), title="test3")
	jim(1:4, 5:9, zeros(4,5), title="test3")
	jim(zeros(4,5), x=1:4, y=5:9, title="test3")
	jim(x=1:4, y=5:9, rand(4,5), title="test4")
	true
end
