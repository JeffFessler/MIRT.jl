#---------------------------------------------------------
# # [ImageGeom](@id ImageGeom)
#---------------------------------------------------------

# This page explains
# [`ImageGeom`](https://github.com/JeffFessler/MIRT.jl/blob/master/src/fbp/image_geom.jl)
# in
# [MIRT.jl](http://github.com/JeffFessler/MIRT.jl).

# ### Setup

# Packages needed here.

using MIRT
using MIRTjim: jim, prompt
using Plots; default(markerstrokecolor=:auto)

# The following line is helpful when running this example.jl file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);

# ### Overview

# When performing tomographic image reconstruction,
# one must specify the geometry of the grid of image pixels.
# (In contrast, in image denoising and image deblurring problems,
# one works with the given discrete image
# and no physical coordinates are needed.)

# The key parameters of a grid of image pixels are
# * the size (dimensions) of the grid, e.g., `128 × 128`,
# * the spacing of the pixels, e.g., `1mm × 1mm`,
# * the offset of the pixels relative to the origin, e.g., `(0,0)`

# The data type `ImageGeom` describes such a geometry,
# for arbitrary dimensions (2D, 3D, etc.).

# There are several ways to construct this structure.
# The default is a `128 × 128` grid
# with pixel size ``\Delta_X = \Delta_Y = 1`` (unitless) and zero offset:

ig = ImageGeom()

# Here is a 3D example with non-cubic voxel size:

ig = ImageGeom( (512,512,128), (1,1,2), (0,0,0) )

# To avoid remembering the order of the arguments,
# named keyword pairs are also supported:

ig = ImageGeom( dims=(512,512,128), deltas=(1,1,2), offsets=(0,0,0) )

# ### Units

# The pixel dimensions `deltas` can (and should!) be values with units.

# Here is an example for a video (2D+time) with 12 frames per second:
using UnitfulRecipes
using Unitful: mm, s

ig = ImageGeom( dims=(640,480,1000), deltas=(1mm,1mm,(1//12)s) )

# ### Properties

# An ImageGeom object has *many* useful properties;
# there are too many to remember, so there is built-in help:

ig.help

# This small 2D example illustrates how some properties are used:

ig = ImageGeom( dims=(12,10), deltas=(1mm,1mm), offsets=(0.5,0.5) )

showgrid = (ig) -> # x,y grid locations of pixel centers
    scatter(ig.xg, ig.yg, label="", xlabel="x", ylabel="y",
        xlims = maximum(abs, ig.x) * 1.2 .* (-1,1),
        xticks=[ig.x[1], zero(eltype(ig.x)), ig.x[end]],
        ylims = maximum(abs, ig.y) * 1.2 .* (-1,1),
        yticks=[ig.y[1], zero(eltype(ig.y)), ig.y[end]],
        aspect_ratio=1, title="offsets $(ig.offsets)")
showgrid(ig)

#-
prompt();


# ### Offsets (unitless translation of grid)

# The default `offsets` are zeros,
# corresponding to symmetric sampling around origin:

ig = ImageGeom( dims=(12,10), deltas=(1mm,1mm) )
p = showgrid(ig)

#-
prompt();

# That default offset is natural for tomography
# when considering finite pixel size:

square = (x,y,Δ) -> plot!(p, label="", color=:black,
    x .+ Δ[1] * ([0,1,1,0,0] .- 0.5),
    y .+ Δ[2] * ([0,0,1,1,0] .- 0.5),
)
showgrid(ig)
square2 = (x,y) -> square(x, y, ig.deltas)
square2.(ig.xg,ig.yg)
plot!(p)

#-
prompt();

# In that default geometry, the center `(0,0)` of the image
# is at a corner of the middle 4 pixels (for even image sizes).
# That default is typical for tomographic imaging (e.g., CT, PET, SPECT).
# One must be careful when using operations like `imrotate` or `fft`.


# ### Mask

# In tomographic image reconstruction, patients are usually more "round"
# than "square" so often we only want to estimate the pixels inside some
# support `mask`: a Bool array indicating which pixels are to be estimated.
# (The rest are constrained to be zero.)
# The `ImageGeom` struct has an entry to store this `mask`.
# The default is `Trues(dims)` which is a "lazy" Bool `AbstractArray`
# from the `FillArrays` package that is conceptually similar to `trues(dims)`
# but requires `O(1)` storage.  So there is essentially no memory penalty
# to storing this entry in the `ImageGeom` for users who do not want
# to think about a `mask`.
# For users who do want a `mask`, fortunately Julia uses a special `BitArray`
# type to store Bool arrays, so the storage is 8× less than using bytes
# in most other languages.

# Often we use a "circle inscribed in the square" as a generic support mask,
# and one of the built-in properties can generate such a circular mask:

mask = ig.circ()
ig = ImageGeom(ig.dims, ig.deltas, ig.offsets, mask)
ig.plot(jim)

# Note that `jim` displays the axes with the units naturally;
# see [MIRTjim.jl](http://github.com/JeffFessler/MIRTjim.jl).

# There is an older interface `image_geom` similar to the function
# of the same name in
# [the Matlab version of MIRT](https://github.com/JeffFessler/mirt)
# provided for backward compatibility,
# but using `ImageGeom` is recommended for Julia work.


# ### AxisArrays

# There is a natural connection between `ImageGeom` and `AxisArrays`.
# Note the automatic labeling by
# [MIRTjim.jim](https://github.com/JeffFessler/MIRTjim.jl).

using AxisArrays
using Unitful: mm
ig = ImageGeom( dims=(60,48), deltas=(1.5mm,1mm) )
za = AxisArray( ig.circ() * 10/mm ; x=ig.x, y=ig.y)
jim(za, "AxisArray example")
