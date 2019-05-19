* Install Julia (1.1 or later) from https://julialang.org/
* Explore the extensive documentation at https://docs.julialang.org/
* Launch Julia and use the `]` key to enter its package manager.
* Add any packages needed for these notebooks using the `add` command.
* For example `add Plots` to add the `Plots` package.
* Other crucial standard packages are `FFTW` `FFTViews`
* You will also need the package `IJulia` to run any of these demo notebooks.
* Add MIRT by typing `add MIRT`
(Automatically loads from https://github.com/JeffFessler/MIRT.jl because MIRT.jl is a registered package.)
* Type `precompile` to have Julia precompile the added packages.
* After you are done adding packages, press the backspace key to return to the REPL prompt.
* Later if you need to add more packages just type the `]` key again at the REPL prompt to enter the package manager.
* Julia is under active development so code is updated frequently.  It is a wise to type `up` (short for `update`) at the package manager prompt every week or so to get automatic updates of any packages you have intalled.
* For some Julia tutorials see
http://web.eecs.umich.edu/~fessler/course/551/julia/tutor/
* For some signal processing demos in Julia see
http://web.eecs.umich.edu/~fessler/course/551/julia/demo/
