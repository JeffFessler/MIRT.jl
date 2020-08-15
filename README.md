# MIRT (Michigan Image Reconstruction Toolbox) in Julia

[![Build Status][action-img]][action-url]
[![Build Status][pkgeval-img]][pkgeval-url]
[![Codecov.io][codecov-img]][codecov-url]
[![Coveralls][coveralls-img]][coveralls-url]
[![docs][docs-stable-img]][docs-stable-url]
[![docs][docs-dev-img]][docs-dev-url]

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
- [doc/start.md](https://github.com/JeffFessler/MIRT.jl/blob/master/doc/start.md)
for conventional Julia
- [doc/start-pro.md](https://github.com/JeffFessler/MIRT.jl/blob/master/doc/start-pro.md)
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

Current version is tested with Julia 1.5.
Older tagged versions should work with older Julia versions.

<!-- URLs -->
[action-img]: https://github.com/JeffFessler/MIRT.jl/workflows/Unit%20test/badge.svg
[action-url]: https://github.com/JeffFessler/MIRT.jl/actions
[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/I/ImageDraw.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html
[travis-img]: https://travis-ci.org/JeffFessler/MIRT.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JeffFessler/MIRT.jl
[codecov-img]: https://codecov.io/github/JeffFessler/MIRT.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/JeffFessler/MIRT.jl?branch=master
[coveralls-img]: https://coveralls.io/repos/JeffFessler/MIRT.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JeffFessler/MIRT.jl?branch=master
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JeffFessler.github.io/MIRT.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JeffFessler.github.io/MIRT.jl/dev
