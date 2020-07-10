# coil_compress.jl

using MIRT: ir_mri_coil_compress
using MIRT: jim, prompt, ir_mri_sensemap_sim, ir_load_brainweb_t1_256
import MIRT: snr2sigma

using LinearAlgebra: svd, norm
using Plots: plot, plot!, scatter, scatter!
using Test: @inferred


xtrue = ir_load_brainweb_t1_256()
(nx, ny) = size(xtrue)
ncoil = 8
smap = ir_mri_sensemap_sim(dims = (nx, ny), dx=192/nx, ncoil=ncoil, rcoil=100)
idata = smap .* xtrue
σ = snr2sigma(50, idata)
idata = idata + σ * randn(ComplexF32, size(idata))
nkeep = 4
#@inferred ir_mri_coil_compress(idata, ncoil = nkeep) # todo: fails
(odata, σ, _) = ir_mri_coil_compress(idata, ncoil = nkeep)
for iz=1:nkeep # normalize for display
	odata[:,:,iz] = odata[:,:,iz] / maximum(abs.(odata[:,:,iz]))
end
p = scatter(1:nkeep, σ[1:nkeep], marker=:circle, label="")
scatter!(p, (nkeep+1):ncoil, σ[nkeep+1:ncoil], marker=:x, label="")
tmp = round((norm(σ[1:nkeep]) / norm(σ))^2*100, digits=2)
plot!(p, title="percent kept $tmp")
plot(jim(smap, "smap"), jim(idata, "idata"), jim(odata, "odata"), p)
