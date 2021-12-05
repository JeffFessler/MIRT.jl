# MIRT (Michigan Image Reconstruction Toolbox) in Julia

[![action status][action-img]][action-url]
[![pkgeval status][pkgeval-img]][pkgeval-url]
[![codecov][codecov-img]][codecov-url]
[![coveralls][coveralls-img]][coveralls-url]
[![docs stable][docs-stable-img]][docs-stable-url]
[![docs dev][docs-dev-img]][docs-dev-url]
[![license][license-img]][license-url]
<img src="deps/mirt-logo.svg" alt="MIRTlogo" width="150"/>

https://github.com/JeffFessler/MIRT.jl

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
for conventional Julia
- [doc/start-pro.md](https://github.com/JeffFessler/MIRT.jl/blob/main/doc/start-pro.md)
for
[JuliaPro](https://juliacomputing.com/products/juliapro.html)

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
http://web.eecs.umich.edu/~fessler/irt/reproduce/


## Compatibility

Current version is tested with "latest" stable version of Julia
(currently 1.6.1).
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


<!-- URLs -->
[action-img]: https://github.com/JeffFessler/MIRT.jl/workflows/Unit%20test/badge.svg
[action-url]: https://github.com/JeffFessler/MIRT.jl/actions
[build-img]: https://github.com/JeffFessler/MIRT.jl/workflows/CI/badge.svg?branch=main
[build-url]: https://github.com/JeffFessler/MIRT.jl/actions?query=workflow%3ACI+branch%3Amain
[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/M/MIRT.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/M/MIRT.html
[codecov-img]: https://codecov.io/github/JeffFessler/MIRT.jl/coverage.svg?branch=main
[codecov-url]: https://codecov.io/github/JeffFessler/MIRT.jl?branch=main
[coveralls-img]: https://coveralls.io/repos/JeffFessler/MIRT.jl/badge.svg?branch=main
[coveralls-url]: https://coveralls.io/github/JeffFessler/MIRT.jl?branch=main
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JeffFessler.github.io/MIRT.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JeffFessler.github.io/MIRT.jl/dev
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[license-url]: LICENSE
