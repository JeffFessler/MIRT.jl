#=
Project:      MIRTdev
File:         mri_geom_compress.jl
Description:  MRI Coil Compression
Author:       Rajas Gupta
Date Created: October 8th 2019
Date Update:  October 8th 2019
=#

export ir_mri_coil_compress

using LinearAlgebra: svd, norm


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
    snr2sigma(db, yb)
Convert SNR in dB to noise σ for complex gaussian noise.
No `sqrt(2)` factors is needed here
because `randn(Complex{Float})`
already accounts for that.
(See `randn` documentation.)
"""
function snr2sigma(db, yb::AbstractArray{<:Complex})
    return 10^(-db/20) * norm(yb) / sqrt(length(yb))
end
