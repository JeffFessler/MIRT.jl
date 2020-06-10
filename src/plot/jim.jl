#=
jim.jl
jiffy image display
2019-02-23 Jeff Fessler, University of Michigan
=#

export jim

using Plots: heatmap, plot!, annotate!
using LaTeXStrings
using MosaicViews: mosaicview
using FFTViews: FFTView
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
	:tickdigit => 1,
	:title => "",
	:xlabel => "",
	:ylabel => "",
	:fft0 => false,
	:yflip => nothing, # defer to minimum value of y - see default below
	:yreverse => nothing, # defer to whether y is non-ascending
	:abswarn => true, # warn when taking abs of complex images?
	])

minfloor = x -> floor(minimum(x), digits=jim_def[:tickdigit])
maxceil = x -> ceil(maximum(x), digits=jim_def[:tickdigit])

"""
    nothing_else(x, y)
return `y` if `x` is nothing, else return `x`
"""
function nothing_else(x, y)
	return x == nothing ? y : x
end


"""
    jim(z, ...)

A jiffy image display of `x` using `heatmap`

in
- `z` image, can be 2D or higher, if higher then it uses `mosaicviews`

option
- `aspect_ratio`; default `:equal`
- `clim`; default `(minimum(z),maximum(z))`
- `color` (colormap, e.g. `:hsv`); default `:grays`
- `ncol` for mosaicview for 3D and higher arrays; default `0` does auto select
- `padval` padding value for mosaic view; default `minimum(z)`
- `line3plot` lines around sub image for 3d mosaic; default `false`
- `line3type` line type around sub image for 3d mosaic; default `(:yellow)`
- `mosaic_npad` # of pixel padding for mosaic view; default `1`
- `fft0` if true use FFTView to display; default `false`
- `title`; default `""`
- `xlabel`; default `""`
- `ylabel`; default `""`
- `yflip`; default `true` if `minimum(y) >= 0`
- `yreverse`; default `true` if `y[1] > y[end]`
- `x` values for x axis; default `1:size(z,1)`
- `y` values for y axis; default `1:size(z,2)`
- `xtick`; default `[minimum(x),maximum(x)]`
- `ytick`; default `[minimum(y),maximum(y)]`

out
- returns plot handle, type `Plots.Plot`

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
		title::AbstractString = jim_def[:title],
		xlabel::AbstractString = jim_def[:xlabel],
		ylabel::AbstractString = jim_def[:ylabel],
		fft0::Bool = jim_def[:fft0],
		x = fft0 ? Int.(-size(z,1)/2:size(z,1)/2-1) : (1:size(z,1)),
		y = fft0 ? Int.(-size(z,2)/2:size(z,2)/2-1) : (1:size(z,2)),
		xtick = (minimum(x) < 0 && maximum(x) > 0) ?
			[minfloor(x),0,maxceil(x)] : [minfloor(x),maxceil(x)],
		ytick = (minimum(y) < 0 && maximum(y) > 0) ?
			[minfloor(y),0,maxceil(y)] : [minfloor(y),maxceil(y)],
		yflip::Bool = nothing_else(jim_def[:yflip], minimum(y) >= 0),
		yreverse::Bool = nothing_else(jim_def[:yreverse], y[1] > y[end]),
		abswarn::Bool = jim_def[:abswarn], # ignored here
	)

	# because some backends require y to be in ascending order
	if yreverse
		y = reverse(y)
		z = reverse(z, dims=2)
	end

	n1,n2,n3 = size(z,1), size(z,2), size(z,3)
	if n3 == 1
		z = z[:,:,1] # revert to 2D
	end
	xy = (x, y)
	if ndims(z) > 2
		n1 += mosaic_npad
		n2 += mosaic_npad
		n3 = size(z,3)
		if ncol == 0
			ncol = Int(floor(sqrt(prod(size(z)[3:end]))))
		end
		z = mosaicview(z ; fillvalue=padval, ncol=ncol, npad=mosaic_npad)
		xy = () # no x,y for mosaic
	elseif fft0
		z = FFTView(z)[x,y]
	end

	if minimum(z) == maximum(z) # uniform image
		x = 1:size(z,1)
		y = 1:size(z,2)
		plot(aspect_ratio=aspect_ratio,
			xlim = [x[1], x[end]],
			ylim = [y[1], y[end]],
			title=title,
			yflip=yflip,
			xlabel=xlabel,
			ylabel=ylabel,
			xtick=xtick,
			ytick=ytick,
		)
		annotate!((sum(x)/length(x), sum(y)/length(y), "Uniform $(z[1])", :red))
	else
		heatmap(xy..., z', transpose=false,
			aspect_ratio=aspect_ratio,
			clim=clim,
			color=color,
			title=title,
			yflip=yflip,
			xlabel=xlabel,
			ylabel=ylabel,
			xtick=xtick,
			ytick=ytick,
		)
	end

	if n3 > 1 && line3plot # lines around each subimage
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
    jim(z, title ; kwargs...)
"""
jim(z::AbstractArray{<:Number}, title::AbstractString ; kwargs...) =
	jim(z ; title=title, kwargs...)


"""
    jim(x, y, z ; kwargs...)
"""
jim(x, y, z ; kwargs...) = jim(z ; x=x, y=y, kwargs...)


"""
    jim(x, y, z, title ; kwargs...)
"""
jim(x, y, z, title::AbstractString ; kwargs...) =
	jim(z ; x=x, y=y, title=title, kwargs...)


"""
    jim()
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
    jim(key, value)
set default value for one of the keys
"""
function jim(key::Symbol, value)
	global jim_def
	!haskey(jim_def, key) && throw(ArgumentError("no key $key"))
	jim_def[key] = value
end


"""
    jim(:test)

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

	jim(ones(4,3), title="test2", xlabel=L"x")
	jim(rand(4,3,5), title=L"test3 x^2_i")
	jim(1:4, 5:9, zeros(4,5), title="test3", ylabel=L"y")
	jim(zeros(4,5), x=1:4, y=5:9, title="test3")
	jim(rand(6,4), fft0=true)
	jim(x=1:4, y=5:9, rand(4,5), title="test4")
	jim(x=-9:9, y=9:-1:-9, (-9:9) * (abs.((9:-1:-9) .- 5) .< 3)', title="rev")
	jim(rand(4,5), color=:hsv)
	jim(ones(3,3)) # uniform
	jim(:abswarn, false)
	jim(complex(rand(4,3)))
	jim(complex(rand(4,3)), "complex")
	jim(:abswarn, true)
	true
end
