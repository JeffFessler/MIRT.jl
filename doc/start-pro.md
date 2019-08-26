## Getting started instructions for using JuliaPro and MIRT

* Install JuliaPro 1.1 (or later) from https://juliacomputing.com/products/juliapro.html
 - Email registration is needed
 - Refer to the Quickstart Guide for prerequisites and for the REPL layout
 - When prompted about the Pkg server, leave blank to use the default.

* Explore the extensive Julia documentation at https://docs.julialang.org/
* Launch JuliaPro and type ```1+2``` at the Julia REPL prompt to verify it works.
* Use the `]` key at the Julia REPL to enter its package manager (pkg).
* Type `?` and press Enter at the Julia `pkg` prompt to peruse the many pkg commands.
* Add needed packages using the `add` command.
* Specifically: `add MIRT Plots IJulia Debugger FFTW FFTViews DSP`
* The package `IJulia` is needed for Jupyter notebooks.
* FYI: `add MIRT` automatically loads from https://github.com/JeffFessler/MIRT.jl because MIRT.jl is a registered package.
* Type `precompile` to have Julia precompile the added packages.
* After you are done adding packages, press the backspace key to return to the REPL prompt.
* Later if you need to add more packages just type the `]` key again at the REPL prompt to enter the package manager.
* Julia is under active development so code is updated frequently.  It is a wise to type `up` (short for `update`) at the package manager prompt every week or so to get automatic updates of any packages you have installed.
* After running `up`, you can type `gc` (for garbage collect) to have Julia remove old versions of packages.
* For some Julia tutorials see
http://web.eecs.umich.edu/~fessler/course/551/julia/tutor/
* For some signal processing demos in Julia see
http://web.eecs.umich.edu/~fessler/course/551/julia/demo/
* Test your `MIRT` installation by displaying a test image:
```
using MIRT: jim, ellipse_im
x = ellipse_im(100)
jim(x, title="test")
```
* This test should produce a grayscale image of the famous
[Shepp-Logan phantom](https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom) in the Plots tab of JuliaPro.
* To learn about the jiffy image display function `jim`, type `?jim` at the REPL.
* JuliaPro uses Juno and Atom and there is lots of online docs about these.
I use the vim-mode-plus key bindings, installed using Atom preferences.
