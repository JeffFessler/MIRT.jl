# MIRT (Michigan Image Reconstruction Toolbox) in Julia

https://github.com/JeffFessler/MIRT.jl
<img src="deps/mirt-logo.svg" alt="MIRTlogo" width="150">

[![docs-stable][docs-stable-img]][docs-stable-url]
[![docs-dev][docs-dev-img]][docs-dev-url]
[![action][action-img]][action-url]
[![Aqua QA][aqua-img]][aqua-url]
[![codecov][codecov-img]][codecov-url]
[![deps][deps-img]][deps-url]
[![license][license-img]][license-url]
[![pkgeval][pkgeval-img]][pkgeval-url]
[![version][ver-img]][ver-url]

This is a collection of tools for
[image reconstruction](https://en.wikipedia.org/wiki/Iterative_reconstruction)
in the open-source
[Julia language](https://julialang.org/).

Currently it contains a limited
collection of the tools from the
[Matlab version of MIRT](https://github.com/JeffFessler/mirt),
but this collection will grow over time.

This software was developed at the
[University of Michigan](https://umich.edu/)
by
[Jeff Fessler](http://web.eecs.umich.edu/~fessler)
and his
[group](http://web.eecs.umich.edu/~fessler/group),
with valuable contributions from the Julia community.


## Getting started

For detailed installation instructions, see:
- [doc/start.md](https://github.com/JeffFessler/MIRT.jl/blob/main/doc/start.md)

This package is registered in the
[`General`](https://github.com/JuliaRegistries/General) registry,
so you can install at the REPL with `] add MIRT`.

For development,
use the `dev` command in Julia's package manager,
or
* `git clone https://github.com/JeffFessler/MIRT.jl`


## Examples

You can test drive some jupyter notebooks in your browser
without installing any local software
by using the free service at
https://mybinder.org/
through the links at the companion demonstration repository
* https://github.com/JeffFessler/mirt-demo


## Reproducible research

This package includes some tools from
https://gitlab.eecs.umich.edu/michigan-fast-optimization

For additional reproducible research code, see
http://web.eecs.umich.edu/~fessler/irt/reproduce


## Compatibility

Tested with Julia ≥ 1.12.
Older tagged versions should work with older Julia versions.


## Related packages

* [JuliaImageRecon](https://github.com/JuliaImageRecon)
  growing collection of image reconstruction packages
* [LinearMapsAA](https://github.com/JeffFessler/LinearMapsAA.jl)
  is central to how imaging system models are used here.


## Deprecations

Early versions of `MIRT.jl`
attempted to house many methods
under one roof.
More recently,
the methods are being isolated
into smaller component packages at
[JuliaImageRecon](https://github.com/JuliaImageRecon).

A similar evolution happened
with
[Images/Images.jl](https://github.com/JuliaImages/Images.jl)
and it is likely
that MIRT will evolve
to be an "umbrella package"
that exports a set of packages
that are useful for image reconstruction.

As of `v0.15`, the following functions are deprecated:
| old | new | see |
| :--- | :--- | :--- |
| `jim` | `MIRTjim.jim` | [MIRTjim.jl](https://github.com/JeffFessler/MIRTjim.jl) |
| `prompt` | `MIRTjim.prompt` | |
| `fld_read` | `FileIO.load` | [AVSfldIO.jl](https://github.com/JeffFessler/AVSfldIO.jl) |
| `fld_write` | `FileIO.save` | [FileIO.jl](https://github.com/JuliaIO/FileIO.jl) |
| `ndgrid` | `LazyGrids.ndgrid` | [LazyGrids.jl](https://github.com/JuliaArrays/LazyGrids.jl) |
| `image_geom` | `ImageGeoms.ImageGeom` | [ImageGeoms.jl](https://github.com/JuliaImageRecon/ImageGeoms.jl) |
| `ellipse_im` | `ImagePhantoms.phantom` | [ImagePhantoms.jl](https://github.com/JuliaImageRecon/ImagePhantoms.jl) |
| `ellipse_sino` | `ImagePhantoms.radon` | [ImagePhantoms.jl](https://github.com/JuliaImageRecon/ImagePhantoms.jl) |
| `mri_objects` | `ImagePhantoms.spectrum` | [ImagePhantoms.jl](https://github.com/JuliaImageRecon/ImagePhantoms.jl) |
| `sino_geom` | `Sinograms.RayGeom` | [Sinograms.jl](https://github.com/JuliaImageRecon/Sinograms.jl) |


<!-- URLs -->
[action-img]: https://github.com/JeffFessler/MIRT.jl/workflows/CI/badge.svg
[action-url]: https://github.com/JeffFessler/MIRT.jl/actions

[aqua-img]: https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[codecov-img]: https://codecov.io/github/JeffFessler/MIRT.jl/coverage.svg
[codecov-url]: https://codecov.io/github/JeffFessler/MIRT.jl

[deps-img]: https://juliahub.com/docs/MIRT/deps.svg
[deps-url]: https://juliahub.com/ui/Packages/MIRT

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JeffFessler.github.io/MIRT.jl/dev
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JeffFessler.github.io/MIRT.jl/stable

[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg
[license-url]: LICENSE

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/M/MIRT.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/M/MIRT.html

[ver-img]: https://juliahub.com/docs/MIRT/version.svg
[ver-url]: https://juliahub.com/ui/Packages/MIRT
