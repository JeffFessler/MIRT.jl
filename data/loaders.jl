# loaders.jl
# utilities for loading in data / images

using FileIO
using Colors

function ir_load_brainweb_t1_256()
	file = dirname(pathof(MIRT)) * "/../data/mri/brainweb_t1_256.gif"
	return Float32.(Gray.(load(file)))' # transpose!
end
