## Converting Matlab to Julia in MIRT.

Many of the functions in MIRT.jl were converted from corresponding
functions in the
[Matlab version of MIRT](https://github.com/JeffFessler/mirt/tree/master)

Studying them can be helpful for planning translations of other Matlab functions.

Here is a list of Julia functions that have direct counter-parts
in Matlab MIRT directories of the same base name:

* `fbp/cuboid_im.jl`
* `fbp/ellipse_im.jl`
* `fbp/ellipse_sino.jl`
* `fbp/ellipsoid_im.jl`
* `fbp/image_geom.jl`
* `fbp/rect_im.jl`
* `fbp/rect_sino.jl`
* `fbp/sino_geom.jl`

Here are some that come from MIRT files in different directories:

* `io/fld-read.jl`	<- `utilities/fld_read.m`
* `io/fld-write.jl`	<- `utilities/fld_write.m`
* `mri/sensemap-sim.jl`	<- `mri/mri_sensemap_sim.m`
* `mri/kspace.jl`	<- `ir_mri_kspace_ga_radial.m`
* `nufft/dtft.jl`	<- `nufft/dtft.m`
* `regularize/Aodwt.jl`	<- `penalty/Godwt1.m`
* `utility/mask.jl`	<- `utility/masker.m` `utility/embed.m`


Matlab does not support named keyword arguments,
so the Matlab version of MIRT uses a custom approach called
[`vararg_pair`](https://github.com/JeffFessler/mirt/blob/master/utilities/vararg_pair.m)
as a work-around.
In Julia you replace any such use
with
[named keyword arguments](https://docs.julialang.org/en/v1/manual/functions/#Keyword-Arguments)
with appropriate default values.

When translating code from Matlab to Julia,
often you will need to locate some functions in Matlab version of MIRT.
To search for those functions,
one way is to go to the repo
https://github.com/JeffFessler/mirt
and then enter the function name in the github search bar.
Github will then return a list like this
https://github.com/JeffFessler/mirt/search?q=vararg_pair
