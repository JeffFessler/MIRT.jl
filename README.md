# MIRT (Michigan Image Reconstruction Toolbox) in Julia

[![Build Status](https://travis-ci.org/JeffFessler/MIRT.jl.svg?branch=master)](https://travis-ci.org/JeffFessler/MIRT.jl) 
[![codecov.io](http://codecov.io/github/JeffFessler/MIRT.jl/coverage.svg?branch=master)](http://codecov.io/github/JeffFessler/MIRT.jl?branch=master) 
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
[group](http://web.eecs.umich.edu/~fessler/group).


## Getting started

For detailed installation instructions, see:
[doc/start.md](https://github.com/JeffFessler/MIRT.jl/blob/master/doc/start.md)

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
