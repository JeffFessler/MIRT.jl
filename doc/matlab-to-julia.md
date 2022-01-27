## Converting Matlab to Julia in MIRT.

Many of the functions in MIRT.jl were converted from corresponding
functions in the
[Matlab version of MIRT](https://github.com/JeffFessler/mirt/tree/main).

Studying them can be helpful
for planning translations of other Matlab functions.
The languages have many similarities,
but also many dramatic differences
(especially multiple dispatch)
so automatic translation
seems impossible
and one must really understand the algorithm
to do a good conversion.

Here is a list of functions
in the
[Julia MIRT `fbp/`][fbp-jl] directory
that have direct counter-parts
in the Matlab
[Matlab MIRT `fbp/`][fbp-mat] directory:

* `fbp/cuboid_im.jl`
* `fbp/ellipse_im.jl`
* `fbp/ellipse_sino.jl`
* `fbp/ellipsoid_im.jl`
* `fbp/image_geom.jl`
* `fbp/rect_im.jl`
* `fbp/rect_sino.jl`
* `fbp/sino_geom.jl`


Here are some pairs
that came from MIRT files in different directories:

| Julia | Matlab |
| :-- | :-- |
| [`AVSfldIO/fld-read.jl`][fld-read-jl] | [`utilities/fld_read.m`][fld-read-mat] |
| [`AVSfldIO/fld-write.jl`][fld-write-jl] | [`utilities/fld_write.m`][fld-write-mat] |
| [`mri/sensemap-sim.jl`][sensemap-sim-jl] | [`mri/mri_sensemap_sim.m`][sensemap-sim-mat] |
| [`mri/kspace.jl`][kspace-jl] | [`mri/ir_mri_kspace_ga_radial.m`][kspace-mat] |
| [`nufft/dtft.jl`][dtft-jl] | [`nufft/dtft.m`][dtft-mat] |
| [`regularize/Aodwt.jl`][odwt-jl] | [`penalty/Godwt1.m`][odwt-mat] |
| [`utility/mask.jl`][mask-jl] | [`utility/masker.m`][masker-mat] & [`utility/embed.m`][embed-mat] |


Matlab does not support named keyword arguments,
so the Matlab version of MIRT uses a custom approach called
[`vararg_pair`][vararg-mat]
as a work-around.
In Julia, replace any such use with
[named keyword arguments](https://docs.julialang.org/en/v1/manual/functions/#Keyword-Arguments)
with appropriate default values.

When translating code from Matlab to Julia,
often you will need to locate some functions in Matlab version of MIRT.
To search for those functions,
one way is to go to the
[mirt repo](https://github.com/JeffFessler/mirt)
and then enter the function name in the github search bar.
Github will then return a list
[like this][vararg-search].


<!-- URLs -->
[vararg-search]: https://github.com/JeffFessler/mirt/search?q=vararg_pair
[vararg-mat]: https://github.com/JeffFessler/mirt/blob/main/utilities/vararg_pair.m
[fbp-jl]: https://github.com/JeffFessler/MIRT.jl/blob/main/src/fbp/
[fbp-mat]: https://github.com/JeffFessler/mirt/blob/main/fbp/
[fld-read-jl]: https://github.com/JeffFessler/AVSfldIO.jl/blob/main/src/fld-read.jl
[fld-read-mat]: https://github.com/JeffFessler/mirt/blob/main/utilities/fld_read.m
[fld-write-jl]: https://github.com/JeffFessler/AVSfldIO.jl/blob/main/src/fld-write.jl
[fld-write-mat]: https://github.com/JeffFessler/mirt/blob/main/utilities/fld_write.m
[dtft-jl]: https://github.com/JeffFessler/MIRT.jl/blob/main/src/nufft/dtft.jl
[dtft-mat]: https://github.com/JeffFessler/mirt/blob/main/nufft/dtft.m
[kspace-jl]: https://github.com/JeffFessler/MIRT.jl/blob/main/src/mri/kspace.jl
[kspace-mat]: https://github.com/JeffFessler/mirt/blob/main/mri/ir_mri_kspace_ga_radial.m
[odwt-jl]: https://github.com/JeffFessler/MIRT.jl/blob/main/src/regularize/Aodwt.jl
[odwt-mat]: https://github.com/JeffFessler/mirt/blob/main/penalty/Godwt1.m
[sensemap-sim-jl]: https://github.com/JeffFessler/MIRT.jl/blob/main/src/mri/sensemap-sim.jl
[sensemap-sim-mat]: https://github.com/JeffFessler/mirt/blob/main/mri/mri_sensemap_sim.m
[mask-jl]: https://github.com/JeffFessler/MIRT.jl/blob/main/src/utility/mask.jl
[masker-mat]: https://github.com/JeffFessler/mirt/blob/main/utilities/masker.m
[embed-mat]: https://github.com/JeffFessler/mirt/blob/main/utilities/embed.m
