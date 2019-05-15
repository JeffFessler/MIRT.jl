# loaders.jl
# utilities for loading in data / images

"""
`data = ir_load_brainweb_t1_256()`

load brainweb T1-weighted MRI slice of size 256x256
"""
function ir_load_brainweb_t1_256()
	file = joinpath(dirname(pathof(MIRT)), "../data/mri/brainweb_t1_256.fld")
	return Float32.(fld_read(file))' # transpose!
end
