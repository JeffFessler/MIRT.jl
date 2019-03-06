# jim.jl

using Plots
using MosaicViews
using FFTViews

# global state variable(s)
jim_state_abswarn = true # warn when taking abs of complex images?


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
* `fft0` if true use FFTView to display; default false
* `title` for heatmap; default `""`
* `xlabel` for heatmap; default `""`
* `ylabel` for heatmap; default `""`
* `yflip` for heatmap; default `true` if minimum(y) >= 0
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
	fft0::Bool = false,
	x = fft0 ? Int.(-size(z,1)/2:size(z,1)/2-1) : (1:size(z,1)),
	y = fft0 ? Int.(-size(z,2)/2:size(z,2)/2-1) : (1:size(z,2)),
	xtick = (minimum(x) < 0 && maximum(x) > 0) ?
		 [minimum(x),0,maximum(x)] : [minimum(x),maximum(x)],
	ytick = (minimum(y) < 0 && maximum(y) > 0) ?
		 [minimum(y),0,maximum(y)] : [minimum(y),maximum(y)],
	yflip = minimum(y) >= 0,
	abswarn::Bool = jim_state_abswarn,
	)

	if !isreal(z)
		if abswarn
			@warn "magnitude"
		end
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
	elseif fft0
		z = FFTView(z)[x,y]
	end

heatmap(x, y, z', transpose=false,
	aspect_ratio=aspect_ratio,
	clim=clim,
	color=color,
	title=title,
	yflip=yflip,
	xlabel=xlabel,
	ylabel=ylabel,
	xtick=xtick,
	ytick=ytick)

end # jim


"""
`jim(z, title; kwargs...)`
"""
function jim(z, title::String; kwargs...)
	return jim(z; title=title, kwargs...)
end


"""
`jim(x, y, z; kwargs...)`
"""
function jim(x, y, z; kwargs...)
	return jim(z; x=x, y=y, kwargs...)
end


"""
`jim(x, y, z, title; kwargs...)`
"""
function jim(x, y, z, title::String; kwargs...)
	return jim(z; x=x, y=y, title=title, kwargs...)
end


# show docstring if user calls it with no arguments
function jim()
	@doc jim
end


"""
`jim(abswarn=false)`
"""
function jim(; abswarn::Bool=jim_state_abswarn)
	global jim_state_abswarn = abswarn
end


"""
`jim(:test)`
"""
function jim(test::Symbol)
	@assert test == :test
	jim(ones(4,3), title="test2")
	jim(ones(4,3,5), title="test3")
	jim(1:4, 5:9, zeros(4,5), title="test3")
	jim(zeros(4,5), x=1:4, y=5:9, title="test3")
	jim(x=1:4, y=5:9, rand(4,5), title="test4")
	true
end
