#=
Project:      MIRTdev
File:         mri_geom_compress.jl
Description:  MRI Coil Compression
Author:       Rajas Gupta
Date Created: October 8th 2019
Date Update:  October 8th 2019
=#

export ir_mri_coil_compress

# using MIRT: jim, prompt, ir_mri_sensemap_sim, ir_load_brainweb_t1_256
using LinearAlgebra: svd, norm
using Plots: plot, plot!, scatter, scatter!
using Test: @inferred


"""
    (odata, σ, Vr) = ir_mri_coil_compress(idata ; ncoil)

MRI coil compression via PCA.
Given multiple MRI surface coil images (idata), use SVD/PCA
to find a smaller number of virtual coil images (odata).

In:
* `idata` `[(N) n_in]`: noisy complex images (2D or 3D) for each coil

Option:
* `ncoil` Desired # of virtual coils (default: 1)

Out:
* `odata` `[(N) ncoil]`: virtual coil images
* `σ`     `[n_in]`: singular values.
* `Vr`    `[n_in, ncoil]`: compression matrix for reducing other data.

todo: currently ignores noise correlations

"""
function ir_mri_coil_compress(idata::AbstractArray{<:Number} ; ncoil::Int = 1)

    idim = size(idata)
    n_in = idim[end]
    idata = reshape(idata, :, n_in)

    (_, σ, V) = svd(idata)

    Vr = V[:, 1:ncoil]
    odata = idata * Vr
    odata = reshape(odata, idim[1:end-1]..., ncoil)
    return (odata, σ, Vr)
end


"""
    ir_mri_coil_compress(:test)
self test with plots
"""
function ir_mri_coil_compress(test::Symbol)
    test != :test && throw(DomainError(test, "Not valid"))
    xtrue = ir_load_brainweb_t1_256()
    (nx, ny) = size(xtrue)
    ncoil = 8
    (smap, _) = ir_mri_sensemap_sim(dims = (nx, ny), dx=192/nx, ncoil=ncoil, rcoil=100)
    idata = smap .* xtrue
    sig = snr2sigma(50, idata)
    idata = idata + sig * (randn(size(idata)) + 1im * randn(size(idata)))
    nkeep = 4
#   @inferred ir_mri_coil_compress(idata, ncoil = nkeep) # todo: fails
    (odata, σ, _) = ir_mri_coil_compress(idata, ncoil = nkeep)
    for iz=1:nkeep # normalize for display
        odata[:,:,iz] = odata[:,:,iz] / maximum(abs.(odata[:,:,iz]))
    end
    p = scatter(1:nkeep, σ[1:nkeep], marker=:circle, label="")
    scatter!(p, (nkeep+1):ncoil, σ[nkeep+1:ncoil], marker=:x, label="")
    tmp = round((norm(σ[1:nkeep]) / norm(σ))^2*100, digits=2)
    plot!(p, title="percent kept $tmp")
    plot(jim(smap, "smap"), jim(idata, "idata"), jim(odata, "odata"), p)
    prompt()
    true
end


"""
    snr2sigma(db, yb)
convert SNR in dB to noise σ for complex gaussian noise
"""
function snr2sigma(db, yb::AbstractArray{<:Complex})
    10^(-db/20) * norm(yb[:]) / sqrt(length(yb)) / sqrt(2)
end

# ir_mri_coil_compress(:test)
