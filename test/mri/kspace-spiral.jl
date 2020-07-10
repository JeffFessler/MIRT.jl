# kspace_spiral.jl

using MIRT: mri_kspace_spiral

using Plots; default(markerstrokecolor=:auto, markersize=1, markershape=:circle)


N = 64

k0, o0, g0 = mri_kspace_spiral() # default 22,-1
for fov in (20,21)
	mri_kspace_spiral( ; fov=fov, Nt=-1, warn_nk=false)
end
k5l,_,g5l = mri_kspace_spiral(nl = 5) # interleaves

plot(xlabel="kx", ylabel="ky", aspect_ratio=1)
p1 = scatter!(k0[:,1], k0[:,2], label = "1-shot spiral")

p2 = plot()
scatter!(g0[:,1], label="gx")
scatter!(g0[:,2], label="gy")

p4 = plot(g5l[:,1,:], label="")
plot!(g5l[:,2,:], label="")

p3 = plot(xlabel="kx", ylabel="ky", aspect_ratio=1, title="5-shot spiral")
for ii=1:5
	scatter!(k5l[:,1,ii], k5l[:,2,ii], label="")
end

plot(p1, p2, p3, p4)
