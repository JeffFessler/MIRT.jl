#=
ellipse_sino.jl
2019-07-12, Helena H
2019-07-13, Jeff Fessler, refactor to use sg.grid
=#

export ellipse_sino

#using MIRT: sino_geom, MIRT_sino_geom, downsample2


"""
`sino = ellipse_sino(sg, ells ; oversample=1, xscale=1, yscale=1)`

Create sinogram projection of one or more ellipses.
Works for any sinogram geometry.

in
- `sg::MIRT_sino_geom,`		sinogram geometry object from `sino_geom()`
- `ells::Matrix`			`[ne 6]` ellipse parameters
							`[centx centy radx rady angle_degrees amplitude]`

options
- `oversample::Int`	oversampling factor for emulating "strips"
						default 1: just 1 ray per detector element
- `xscale::Int`		use -1 to flip in x (not recommended); default 1
- `yscale::Int`		use -1 to flip in y (not recommended); default 1

out
- `sino`		`[nb na]` sinogram

To get the sample locations, use `(pos,ang) = sg.grid`
"""
function ellipse_sino(sg::MIRT_sino_geom, ells::AbstractMatrix{<:Real} ;
		oversample::Int = 1, kwargs...)

	sg = sg.over(oversample)
	(rg, ϕg) = sg.grid
	sino = ellipse_sino(rg, ϕg, ells; kwargs...)
	if oversample > 1
		sino = downsample2(sino, [oversample, 1])
	end
	return sino
end


"""
`sino = ellipse_sino(rg::AbstractArray{<:Real}, ϕg::AbstractArray{<:Real},
		ells ; oversample=1, xscale=1, yscale=1)`

Create sinogram projection of one or more ellipses,
for arbitrary radial/angular sampling grid locations `(rg, ϕg)`

in
- `rg::AbstractArray{<:Real}`	radial sampling locations
- `ϕg::AbstractArray{<:Real}`	angular sampling locations (radians)
- `ells::Matrix`			`[ne 6]` ellipse parameters
							`[centx centy radx rady angle_degrees amplitude]`

options
- `xscale::Int`		use -1 to flip in x (not recommended); default 1
- `yscale::Int`		use -1 to flip in y (not recommended); default 1

out
- `sino::AbstractArray{Float32}` same size as `rg` and `ϕg`
"""
function ellipse_sino(rg::AbstractArray{<:Real}, ϕg::AbstractArray{<:Real},
		ells::AbstractMatrix{<:Real} ;
		xscale::Int = 1,
		yscale::Int = 1,
	)

	size(ells,2) != 6 && throw("6 parameters per ellipse")
	size(rg) != size(ϕg) && throw("rg and ϕg size mismatch")

	sino = zeros(Float32, size(rg))

	cangs = cos.(ϕg)
	sangs = sin.(ϕg)

	#loop over ellipses
	#ticker reset
	ne = size(ells, 1)
	for ie in 1:ne
		#ticker(mfilename, ie, ne)
		ell = ells[ie, :]

		cx = ell[1] * xscale
		cy = ell[2] * yscale
		rx = ell[3]
		ry = ell[4]
		(rx <= 0) || (ry <= 0) && throw("need positive radii")
		eang = deg2rad(ell[5])
		val = ell[6]

		if yscale == -1
			eang = -eang
		end
		if xscale == -1
			eang = pi - eang
		end
		scale = 2 * val * rx * ry

		# square of projected radius:
		rp2 = @. (rx * (cangs * cos(eang) + sangs * sin(eang)))^2 +
				(ry * (sangs * cos(eang) - cangs * sin(eang)))^2
		sp = cx * cangs + cy * sangs # radial shift
		dis2 = (rg - sp).^2 # square of distances from center
		@. sino += scale / rp2 * sqrt(max(rp2 - dis2, 0))
	end

	return sino
end
