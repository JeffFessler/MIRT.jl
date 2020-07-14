# Usage:
#     julia benchmark/run_benchmarks.jl
using BenchmarkTools
using MIRT
using FFTW: fft
using Random
using LinearMapsAA: LinearMapAA

include("benchmark_utils.jl")

on_CI = haskey(ENV, "GITHUB_ACTIONS")

const SUITE = BenchmarkGroup()

for benchmark_file in [
    "mri.jl",
    "nufft.jl",
    "regularize.jl",
    "utility.jl",
]
    include(benchmark_file)
end
