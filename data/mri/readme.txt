brainweb_t1_256.gif
came from brainwed database then in Julia:
	using FileIO
	x = load("brainweb_t1.gif")
	y = x[2:end-1,2:end-1]
	save("brainweb_t1_256.gif", y)

to load it:
x = Float32.(Gray.(load("brainweb_t1_256.gif")))
or simply:
x = ir_load_brainweb_t1_256()
