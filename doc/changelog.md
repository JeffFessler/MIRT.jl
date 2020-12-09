In reverse chronological order

2020-12-08
towards 0.13.2 after NFFT.jl update requires `NFFTPlan` to `plan_nfft`.
todo: expose cuda option!

2020-07-23
Breaking change: SinoGeom

Many changes to improve codecov, including more complete `fld_read fld_write`

2020-07-21
Breaking change: ImageGeom

0.12.0 2020-07-14 for real this time

Use tuples in `downsample2` and related functions.
This is a slightly breaking change.

0.12.0 2020-07-13 - false start

0.11.0 2020-07-13 (because it would not let me tag 0.12 here)

2020-07-10
Breaking change for `ir_mri_sensemap_sim` to isolate plots

2020-07-09
Mildly breaking change: isolate almost all tests to `test/`

0.11.0 2020-07-01 - forgot to tag it

Breaking change: make `LinearMapAO` the default
for `Afft, Anufft, diffl_map, Aodwt`
Also use `dot in `ncg` and `ogm_ls` algorithms
to support array variables.

0.10.0 2020-06-30

2020-06-30
offer LinearMapAO options for `Afft, Anufft, diffl_map`

2020-06-29
embed! getindex!

2020-06-22
In place left finite differences: `diffl!` `diffl_adj!` `diffl_map`

0.9.6 2020-06-10

n-d version of `diff_map`

2020-06-05 isolate `@btime` in time/

0.9.5 2020-04-26

cleanup `:hsv`

0.9.4 2020-04-24

2020-03-21
`mri_trajectory`

0.9.2 2020-02-29

2020-02-17
`reverser`

2020-02-09
`exp_mult`

2020-02-05
`eql_root` `rmsd100`

0.9.1 2020-01-20

`display` to `show`, updated `LinearMapsAA`

0.9.0 2019-12-26 for Julia 1.3

2019-11-22
`mri_image_geom`

0.8.0 2019-11-05

2019-11-07
to Julia 1.2 due to file dependency issues

2019-11-04
`interp1`

2019-08-26
use `LinearMapAA`

2019-07-28
use isolated `MIRTio`

0.7.0 2019-07-15

2019-07-12
add `ellipsoid_im` `ellipse_sino` `rect_sino`

2019-07-07
add `Afft`

2019-07-02
add `cuboid_im` `sino_geom`

2019-06-29
add travis.yml, improve codecov
add `rect_im` `caller_name`
isolate `ndgrid` `rotate2d`

0.6.0 2019-06-23

add `prompt` `ir_mri_sense_map_sim`

2019-06-22
major update to `image_geom`

0.5.0 2019-06-16

add `block_lm`

0.4.0 2019-06-12

2019-06-12
add `dtft` `nufft`

2019-06-10
add `poweriter`

0.3.0 2019-05-23

add `loadpfile` for reading GE MRI scan Pfile (raw k-space data)
fix minor path problem

0.2.0 2019-05-14

add `fld_read` and `fld_write`
