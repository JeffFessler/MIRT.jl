## Converting Matlab to Julia in MIRT.

Many of the functions in MIRT.jl were converted from corresponding
functions in the
[Matlab version of MIRT](https://github.com/JeffFessler/mirt/tree/master)

Studying this can be helpful for planning translations of other Matlab functions.

Here is a list of Julia functions that have direct counter-parts
in Matlab MIRT directories of the same name:

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
