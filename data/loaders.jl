# loaders.jl
# utilities for loading in data / images

using FileIO: load

"""
    data = ir_load_brainweb_t1_256()

Load brainweb T1-weighted MRI slice of size `256 × 256`
"""
function ir_load_brainweb_t1_256()
	dirmirt = normpath(dirname(pathof(MIRT)), "..")
	file = joinpath(dirmirt, "data", "mri", "brainweb_t1_256.fld")
	return Float32.(load(file))' # transpose!
end
