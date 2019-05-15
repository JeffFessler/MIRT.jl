brainweb_t1_256.fld
came from brainweb database http://www.bic.mni.mcgill.ca/brainweb/

to load it simply say
	x = ir_load_brainweb_t1_256()

How formed in Julia (new way 2019-05-14):

using FileIO
x = load("brainweb_t1.gif")
y = x[2:end-1,2:end-1]
save("brainweb_t1_256.gif", y)

# x = Float32.(Gray.(load("brainweb_t1_256.gif"))) # old way

t1 = load("brainweb_t1_256.gif");
t2 = map(f -> f.r, t1);
t3 = t2 * 255;
t4 = UInt8.(round.(t3));
extrema(t4)
extrema(t3) # 0 241
extrema(t4-t3) # 1e-5
fld_write("brainweb_t1_256.fld", t4)
