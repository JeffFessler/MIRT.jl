## Getting started instructions for using MIRT with Julia and Juno/Atom

* Follow the instructions at
http://docs.junolab.org/latest/man/installation/
to set up the [Juno IDE](https://junolab.org/) in
  the [Atom editor](https://atom.io), including:

  - Install [Julia](https://julialang.org) (1.4 or later recommended).
  - Power users: set a shell alias for "julia" to the julia executable.
On my Mac my path is
`~/../freeware/Julia-1.4.app/Contents/Resources/julia/bin/julia`
  - Install the [Atom editor](https://atom.io).
  - Use Atom to install the `uber-juno` package. This may take a couple minutes.
  - Click "Yes" for the question about "Juno-specific panes on startup", unless you are an Atom guru.
  - You should see a "REPL" box at the bottom of Atom.  Hit enter to start Julia.
  - If Julia does not start then follow the [instructions](http://docs.junolab.org/latest/man/installation/#Note-1) to set the Julia path

* Explore the extensive Julia documentation at https://docs.julialang.org/
* Type ```1+2``` at the Julia REPL prompt to verify it works.
* Use the `]` key at the Julia REPL to enter its package manager (pkg).
* Type `?` and press Enter at the Julia `pkg` prompt to peruse the many pkg commands.
* Add needed packages using the `add` command at the `pkg` prompt:
* `add MIRT Plots IJulia FFTW FFTViews Arpack Debugger`
  - The package `IJulia` is needed for Jupyter notebooks.
  - FYI: `add MIRT` automatically loads from https://github.com/JeffFessler/MIRT.jl because MIRT.jl is a registered package.
  - There are numerous packages available, see https://pkg.julialang.org/docs/
* Type `precompile` to have Julia precompile the added packages.
(This will save time later.)
* After you are done adding packages, press the backspace key to return to the REPL prompt.
* Later if you need to add more packages just type the `]` key again at the REPL prompt to enter the package manager.
  - Julia is under active development so code is updated frequently.  It is wise to type `up` (short for `update`) at the package manager prompt every week or so to get automatic updates of any packages you have installed.
  - After running `up`, you can type `gc` (for garbage collect) to have Julia remove old versions of packages.
* For some Julia tutorials see
http://web.eecs.umich.edu/~fessler/course/551/julia/tutor/
* For some signal processing demos in Julia see
http://web.eecs.umich.edu/~fessler/course/551/julia/demo/
* Test your `MIRT` installation by displaying a test image by entering the following code either in the REPL or in the Atom text editor:
```
using MIRT: jim, ellipse_im
x = ellipse_im(100)
jim(x, title="test")
```
* This test should produce a grayscale image of the famous
[Shepp-Logan phantom](https://en.wikipedia.org/wiki/Shepp%E2%80%93Logan_phantom) in the Plots tab of Atom.
* To learn about the jiffy image display function `jim`, type `?jim` at the REPL.
* Juno and Atom have lots of online documention.
I use the `vim-mode-plus` key bindings, installed using Atom preferences.

* To start a Jupyter notebook for Julia, type at the REPL:
`using IJulia; notebook(detached=true)`

* If that fails, then you might need to separately install [Jupyter Notebook](https://jupyter.readthedocs.io/en/latest/install.html)
(there are many ways to do it) and then restart Atom and try again.

* To use the debugger in Atom/Juno, see
[this debugging tutorial](http://docs.junolab.org/latest/man/debugging)

* The default setting in Julia seems to be to use just one thread.
Even your laptop probably has multiple cores and using them all
will help the code run faster.
To use all cores,
set the
[JULIA_NUM_THREADS environment variable](https://docs.julialang.org/en/latest/manual/environment-variables/#JULIA_NUM_THREADS-1)
before starting Julia.
Here are
[Windows instructions](./thread-pc.md).
