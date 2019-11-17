#=
jim.jl
jiffy image display
2019-02-23 Jeff Fessler, University of Michigan
=#

export jim

using Plots: heatmap, ColorGradient
using LaTeXStrings
using MosaicViews: mosaicview
using FFTViews: FFTView
using Colors: HSV
using Test: @test

# global default key/values
jim_def = Dict([
	:aspect_ratio => :equal,
	:clim => nothing,
	:color => :grays,
	:line3plot => true, # lines around sub image for 3d mosaic?
	:line3type => (:yellow),
	:ncol => 0,
	:padval => nothing,
	:mosaic_npad => 1,
	:title => "",
	:xlabel => "",
	:ylabel => "",
	:fft0 => false,
	:yflip => nothing, # defer to minimum value of y - see default below
	:abswarn => true, # warn when taking abs of complex images?
	])


"""
`nothing_else(x, y)`
return `y` if `x` is nothing, else return `x`
"""
function nothing_else(x, y)
	return x == nothing ? y : x
end


"""
`jim(z, ...)`

jiffy image display of `x` using `heatmap`

in
- `z` image, can be 2D or higher, if higher then it uses `mosaicviews`

option
- `aspect_ratio` for heatmap; default `:equal`
- `clim` for heatmap; default `(minimum(z),maximum(z))`
- `color` colormap; default `:grays`
- `ncol` for mosaicview for 3D and higher arrays; default `0` does auto select
- `padval` padding value for mosaic view; default `minimum(z)`
- `line3plot` lines around sub image for 3d mosaic; default `false`
- `line3type` line type around sub image for 3d mosaic; default `(:yellow)`
- `mosaic_npad` # of pixel padding for mosaic view; default `1`
- `fft0` if true use FFTView to display; default `false`
- `title` for heatmap; default `""`
- `xlabel` for heatmap; default `""`
- `ylabel` for heatmap; default `""`
- `yflip` for heatmap; default `true` if `minimum(y) >= 0`
- `x` for x axis; default `1:size(z,1)`
- `y` for y axis; default `1:size(z,2)`
- `xtick` for heatmap; default `[minimum(x),maximum(x)]`
- `ytick` for heatmap; default `[minimum(y),maximum(y)]`

out
- returns plot handle

2019-02-23 Jeff Fessler, University of Michigan
"""
function jim(z::AbstractArray{<:Real} ;
		aspect_ratio = jim_def[:aspect_ratio],
		clim = nothing_else(jim_def[:clim], (minimum(z), maximum(z))),
		color = jim_def[:color],
		line3plot = jim_def[:line3plot],
		line3type = jim_def[:line3type],
		ncol::Int = jim_def[:ncol],
		padval = nothing_else(jim_def[:padval], minimum(z)),
		mosaic_npad::Int = jim_def[:mosaic_npad],
		title::String = jim_def[:title],
		xlabel::String = jim_def[:xlabel],
		ylabel::String = jim_def[:ylabel],
		fft0::Bool = jim_def[:fft0],
		x = fft0 ? Int.(-size(z,1)/2:size(z,1)/2-1) : (1:size(z,1)),
		y = fft0 ? Int.(-size(z,2)/2:size(z,2)/2-1) : (1:size(z,2)),
		xtick = (minimum(x) < 0 && maximum(x) > 0) ?
			[minimum(x),0,maximum(x)] : [minimum(x),maximum(x)],
		ytick = (minimum(y) < 0 && maximum(y) > 0) ?
			[minimum(y),0,maximum(y)] : [minimum(y),maximum(y)],
		yflip::Bool = nothing_else(jim_def[:yflip], minimum(y) >= 0),
		abswarn::Bool = jim_def[:abswarn], # ignored here
	)

	n1,n2,n3 = size(z,1), size(z,2), size(z,3)
	xy = (x, y)
	if ndims(z) > 2
		n1 += mosaic_npad
		n2 += mosaic_npad
		n3 = size(z,3)
		if ncol == 0
			ncol = Int(floor(sqrt(prod(size(z)[3:end]))))
		end
		z = mosaicview(z, padval, ncol=ncol, npad=mosaic_npad)
		xy = () # no x,y for mosaic
	elseif fft0
		z = FFTView(z)[x,y]
	end

	# attempt an HSV colormap for phase images
	if color === :hsv
		color = ColorGradient([HSV(h,1,1) for h=LinRange(0,350,351)])
	end

	heatmap(xy..., z', transpose=false,
		aspect_ratio=aspect_ratio,
		clim=clim,
		color=color,
		title=title,
		yflip=yflip,
		xlabel=xlabel,
		ylabel=ylabel,
		xtick=xtick,
		ytick=ytick)

	if line3plot # lines around each subimage
		m1 = (1+size(z,1)) / n1 # add one because of mosaicview non-edge
		m2 = (1+size(z,2)) / n2
		plot_box! = (ox, oy) -> plot!(
			ox .+ [0,1,1,0,0]*n1 .+ 0.0,
			oy .+ [0,0,1,1,0]*n2 .+ 0.0,
			line=jim_def[:line3type], label="")
		for ii=0:n3-1
			i1 = mod(ii, m1)
			i2 = floor(ii / m1)
			plot_box!(i1*n1, i2*n2)
		end
	end
	plot!()

end # jim


# handle case of complex image
function jim(z::AbstractArray{<:Number} ;
			abswarn::Bool = jim_def[:abswarn],
			kwargs...
		)

	if !(eltype(z) <: Real)
		abswarn && (@warn "magnitude at $(caller_name())")
		z = abs.(z)
	end
	jim(z ; kwargs...)
end


"""
`jim(z, title ; kwargs...)`
"""
function jim(z::AbstractArray{<:Number}, title::Union{String,LaTeXString}
		; kwargs...)
	return jim(z ; title=title, kwargs...)
end


"""
`jim(x, y, z ; kwargs...)`
"""
function jim(x, y, z ; kwargs...)
	return jim(z ; x=x, y=y, kwargs...)
end


"""
`jim(x, y, z, title ; kwargs...)`
"""
function jim(x, y, z, title::Union{String,LaTeXString} ; kwargs...)
	return jim(z ; x=x, y=y, title=title, kwargs...)
end


"""
`jim()`
return docstring if user calls `jim()` with no arguments
"""
function jim()
	@doc jim
end


#"""
#`jim(abswarn=false)`
#"""
#function jim( ; abswarn::Bool=jim_state_abswarn)
#	global jim_state_abswarn = abswarn
#end

"""
`jim(key, value)` set default value for one of the keys
"""
function jim(key::Symbol, value)
	global jim_def
	!haskey(jim_def, key) && throw(ArgumentError("no key $key"))
	jim_def[key] = value
end


"""
`jim(:test)`

`jim(:keys)` return default keys

`jim(:defs)` return `Dict` of default keys / vals

`jim(:key)` return `Dict[key]` if possible
"""
function jim(test::Symbol)
	global jim_def
	if test === :keys
		return keys(jim_def)
	end
	if test === :defs
		return jim_def
	end
	if haskey(jim_def, test)
		return jim_def[test]
	end
	test != :test && throw("symbol $test")
	jim()
	jim(:keys)
	jim(:clim)
	@test typeof(jim(:defs)) <: Dict

	jim(ones(4,3), title="test2")
	jim(ones(4,3,5), title="test3")
	jim(1:4, 5:9, zeros(4,5), title="test3")
	jim(zeros(4,5), x=1:4, y=5:9, title="test3")
	jim(zeros(4,6), fft0=true)
	jim(x=1:4, y=5:9, rand(4,5), title="test4")
	jim(rand(4,5), color=:hsv)
	jim(:abswarn, false)
	jim(complex(rand(4,3)))
	jim(complex(rand(4,3)), "complex")
	jim(:abswarn, true)
	true
end
