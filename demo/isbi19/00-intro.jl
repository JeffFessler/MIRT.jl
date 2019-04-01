# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.5
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

# ### MRI reconstruction introductory notebook for ISBI 2019 tutorial
#
# This jupyter notebook is in a directory with other notebooks
# that demonstrate MRI reconstruction using the Julia language
# and the Michigan Image Reconstruction Toolbox (MIRT).
#
# 2019-03-14 Jeff Fessler, University of Michigan

# ### Initial steps
#
# See getting started instructions here:
# https://github.com/JeffFessler/MIRT.jl/blob/master/doc/start.md

# ### Test your installation
#
# The next cell loads the packages needed for this initial test.
# It might take a little bit of time to run the first time you try it because Julia is essentially a compiled language under the hood, even though it feels like an interactive language to the user, so it will be compiling things behind the scenes.

# load all packages needed for this demo 
using MIRT # https://github.com/JeffFessler/MIRT.jl
using FFTW

# ### View the Shepp-Logan image and its spectrum
# The `ellipse_im` function in MIRT can generate phantom images consisting of ellipses such as the Shepp-Logan image.
#
# The `jim` function in MIRT is a jiffy image display routine.

image = ellipse_im(256, oversample=2)
jim(image, "Shepp-Logan Phantom")

# calculate k-space data using fft
kspace = fft(image) # this is a 2D FFT because image is a 2D array
jim(kspace, fft0=true) # show k-space with 0 at center, ala matlab's "fftshift"

# it's hard to see much on a linear scale, so let's use a log scale
logger = (x;min=-6) -> log10.(max.(abs.(x) / maximum(abs.(x)), (10.)^min))
jim(logger(kspace), "k-space data on log scale", fft0=true)

# ### Getting help
# In addition to the documentation linked above,
# to learn about any command in Julia just type `?` followed by the command name
# and Julia will return the documentation string
# that goes with that command.
# Here is an example.
# This example illustrates the very useful "multiple dispatch" feature of Julia;
# essentially there are multiple versions of this function
# that you can call in several different ways
# with different argument combinations,
# and with optimal keyword arguments.

?jim

# ### End of introduction!
#
# If you got this far you are ready for more interesting demos next.
