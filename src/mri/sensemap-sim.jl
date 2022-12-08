#=
sensemap-sim.jl
Simulate coil sensitivity maps

Matlab notes:
- 2005-6-20, Jeff Fessler and Amanda Funai, University of Michigan
- 2014-08-19 JF more testing, verifying phase is correct
- 2014-09-09 modified for 3D by Mai Le
- 2016-05-03 JF fixes

2019-06-18, Jeff Fessler, University of Michigan
=#

export ir_mri_sensemap_sim

using SpecialFunctions: ellipk, ellipe

#using MIRT: ndgrid


"""
    smap = ir_mri_sensemap_sim(...)

Simulate 2D or 3D sensitivity maps for sensitivity-encoded MRI
based on
[grivich:00:tmf](http://doi.org/10.1119/1.19461).

This code makes maps for multiple coils,
but does not model coupling between coils,
so most likely it is an approximation at best.

# option
- `dims::Dims` image size; default (64, 64)
- `dx::Real` pixel/voxel dimension; default: 3
- `dy::Real` pixel/voxel dimension; default: `dx`
- `dz::Real` ""
- `ncoil::Int` # of coils total; default 4
- `nring::Int` # of rings of coils; default 1
- `rcoil::Real` coil radius; default `dx * nx / 2 * 0.50`
- `dz_coil` ring spacing in z; default `nz*dz/nring`
   * (3D geometry is a cylinder)
- `coil_distance::Real` distance of coil center from isocenter
   * for central ring of coils as a multiple of `FOVx`,
   * where `FOVx=nx*dx`; default 1.2
- `orbit::Real` default 360 [degrees]
- `orbit_start::AbstractVector{<:Real} = fill(0, nring)` [degrees]
- `scale::Symbol`
  + `:none` (default)
  + `ssos_center` make SSoS of center = 1

# out
- `smap [dims ncoil]` simulated sensitivity maps (complex!)

All length parameters must have same units (e.g., mm or cm)
"""
function ir_mri_sensemap_sim( ; kwargs...)
    return ir_mri_sensemap_sim(
        Vector{Tuple{Int,Int}}(undef, 0) ; # ir_ic_pair
        kwargs...,
    )[1]
end


"""
    (smap,info) = ir_mri_sensemap_sim( :all ; kwargs)
Like `ir_mri_sensemap_sim` but also returns `info` with data for all coils,
mainly for testing and plotting.
"""
function ir_mri_sensemap_sim(all::Symbol ;
    nring::Int = 1, ncoil::Int = 4, kwargs...,
)
    all == :all || throw("bad symbol argument $all")
    return ir_mri_sensemap_sim(
        vec(collect(Iterators.product(1:nring,1:ncoil))) ; # ir_ic_pair
        nring=nring, ncoil=ncoil, kwargs...,
    )
end


# handle (usual) scalar orbit_start
#=
(no: causes Method overwritten warning)
function ir_mri_sensemap_sim( ir_ic_pair::Vector{Tuple{Int,Int}} ;
    nring::Int = 1, orbit_start::Real = 0, kwargs...,
)
    return ir_mri_sensemap_sim( ir_ic_pair ; nring = nring,
        orbit_start = repeat([orbit_start], nring),
    )
end
=#


"""
    (smap,info) = ir_mri_sensemap_sim( ir_ic_pair ; kwargs)
Like `ir_mri_sensemap_sim` but also returns `info` with data for specific coils
where `ir_ic_pair::Vector{Tuple{Int,Int}}`.
(Usually used internally only.)
- `info::NamedTuple` geometry information for plots
"""
function ir_mri_sensemap_sim(
    ir_ic_pair::Vector{Tuple{Int,Int}} ; # usually length 0
    dims::Dims = (64,64), # 2D default
    dx::Real = 3,
    dy::Real = dx,
    dz::Real = 3,
    ncoil::Int = 4,
    nring::Int = 1,
    rcoil::Real = dx * dims[1] / 2 * 0.5,
    dz_coil::Real = ((length(dims) == 3) ? dims[3] : 1) * dz / nring,
    coil_distance::Real = 1.2, # multiplies fov/2
    orbit::Real = 360,
    orbit_start::AbstractVector{<:Real} = fill(0,nring),
    scale::Symbol = :none, # or :ssos_center
)

    (2 <= length(dims) <= 3) || throw("2D or 3D only")
    nx = dims[1]
    ny = dims[2]
    nz = (length(dims) == 3) ? dims[3] : 0

    coils_per_ring = round(Int, ncoil / nring)
    ncoil != nring * coils_per_ring && throw("nring must be divisor of ncoil")

    (smap, info) = ir_mri_sensemap_sim_do(
            nx, ny, nz,
            dx, dy, dz,
            ncoil, coils_per_ring, rcoil, dz_coil,
            orbit, orbit_start, coil_distance, ir_ic_pair)

    scale_center = (nz == 0) ?
        1 / sqrt(sum(abs2.(smap[nx÷2,ny÷2,:]))) :
        1 / sqrt(sum(abs2.(smap[nx÷2,ny÷2,nz÷2,:])))

    if scale === :ssos_center
        smap *= Float32(scale_center)
    elseif scale != :none
        throw("scale $scale")
    end

    return (smap, info)
end


"""
    (smap, info) = ir_mri_sensemap_sim_do()
"""
function ir_mri_sensemap_sim_do(
    nx, ny, nz,
    dx, dy, dz, ncoil, ncoilpr, rcoil, dz_coil,
    orbit, orbit_start, coil_distance, ir_ic_pair,
)

    T = Float32
    nring = ncoil ÷ ncoilpr
    rlist = T(rcoil) * ones(T, ncoilpr, nring) # coil radii

    plist = zeros(T, ncoilpr, nring, 3) # position of coil center [x y z]
    nlist = zeros(T, ncoilpr, nring, 3) # normal vector (inward) from coil center
    olist = zeros(T, ncoilpr, nring, 3) # unit vector orthogonal to normal vector in x-y
    ulist = zeros(T, ncoilpr, nring, 3) # upward vector

    length(orbit_start) != nring && throw("bad orbit_start length")

    # cylindrical coil configuration, like abdominal coils
    alist = deg2rad.(orbit) * (0:(ncoilpr-1)) / ncoilpr # coil angles [radians]
    z_ring = ((1:nring) .- (nring+1)/2) * dz_coil
    for ir in 1:nring, ic in 1:ncoilpr
        phi = alist[ic] + deg2rad(orbit_start[ir])
        Rad = max(nx/2 * dx, ny/2 * dy) * coil_distance
        plist[ic,ir,:] = [Rad * [cos(phi), sin(phi)]; z_ring[ir]]
        nlist[ic,ir,:] = -[cos(phi), sin(phi), 0*z_ring[ir]] # cylinder
        olist[ic,ir,:] = [-sin(phi), cos(phi), 0]
        ulist[ic,ir,:] .= [0, 0, 1]
    end

    # object coordinates
    x = T.(((1:nx) .- (nx+1)/2) * dx)
    y = T.(((1:ny) .- (ny+1)/2) * dy)
    z = (nz > 0) ? T.(((1:nz) .- (nz+1)/2) * dz) : [0]
    (xx, yy, zz) = ndgrid(x,y,z)

    data_per_coil = Array{Any}(undef, nring, ncoilpr)

    smap = zeros(Complex{T}, nx, ny, max(nz,1), ncoilpr, nring)
    for ir in 1:nring, ic in 1:ncoilpr
        # rotate coordinates to correspond to coil orientation
        zr =    (xx .- plist[ic,ir,1]) .* nlist[ic,ir,1] +
                (yy .- plist[ic,ir,2]) .* nlist[ic,ir,2] +
                (zz .- plist[ic,ir,3]) .* nlist[ic,ir,3]
        xr = xx .* nlist[ic,ir,2] - yy .* nlist[ic,ir,1]
        yr = zz .- plist[ic,ir,3] # translate along object z axis

        # compute sensitivity vectors in coil coordinates
        tmp = ir_mri_smap1.(xr, yr, zr, rlist[ic,ir])
        sx = [p[1] for p in tmp]
        sy = [p[2] for p in tmp]
        sz = [p[3] for p in tmp]

        # coil response depends on transverse magnetization only?
        # todo: unsure if this should depend on sy and ulist in 3D
            bx = sz * nlist[ic,ir,1] + sx * olist[ic,ir,1]
            by = sz * nlist[ic,ir,2] + sx * olist[ic,ir,2]
        #    bz = sz * nlist[ic,ir,3] + sx * olist[ic,ir,3]
        smap[:,:,:,ic,ir] = complex.(bx, by)

        if (ir,ic) ∈ ir_ic_pair # save data for plotting
            data_per_coil[ir, ic] =
                (xr=xr, yr=yr, zr=zr, sx=sx, sy=sy, sz=sz, bx=bx, by=by)
        end
    end

    smap *= rlist[1] / T(2π) # trick: scale so maximum is near unity
    smap = nz == 0 ?
        reshape(smap, nx, ny, ncoil) :
        reshape(smap, nx, ny, nz, ncoil)

    info = (x=x, y=y, z=z, dx=dx, dy=dy, dz=dz,
        nlist=nlist, plist=plist, rlist=rlist, olist=olist, ulist=ulist,
        nring=nring, ncoilpr=ncoilpr, rcoil=rcoil, data=data_per_coil)
    return smap, info
end


"""
    ir_mri_smap_r(r, z)
Function for testing near 0.
"""
function ir_mri_smap_r(r, z)
    M = 4 * r / ((1 + r)^2 + z^2) # = k^2, see ellipke
    (K,E) = (ellipk.(M), ellipe.(M))
    return 2 * z / r * ((1 + r)^2 + z^2)^(-0.5) *
        ((1 + r^2 + z^2) / ((1 - r)^2 + z^2) * E - K)
end


"""
    ir_mri_smap1()

Based on grivich:00:tmf
for a circular coil in "x-y plane" of radius "a"

Note that coil x-y plane is not same as object x-y plane!

Returns `(i,j,k)` components of ``B`` vector for each `(x,y,z)` location.
"""
function ir_mri_smap1(x, y, z, a)
    x = x / a # normalized units
    y = y / a
    z = z / a
    r = sqrt(x^2 + y^2)
    M = 4 * r / ((1 + r)^2 + z^2) # = k^2, see ellipke
    (K,E) = (ellipk.(M), ellipe.(M))

    # the following is B_z in eqn (18) in grivich:00:tmf
    # and same as eqn [10] in wang:00:dop to within constant scale factor
    smap_z = 2 * ((1 + r)^2 + z^2)^(-0.5) *
        (K + (1 - r^2 - z^2) / ((1 - r)^2 + z^2) * E)
    smap_z /= a

    # the following is B_r in eqn (17) in grivich:00:tmf
    smap_r = 2 * z / r * ((1+r)^2 + z^2)^(-0.5) *
        ((1 + r^2 + z^2) / ((1-r)^2 + z^2) * E - K)

    if abs(r) < 1e-6
        smap_r = 3 * pi * z / ((1 + z^2).^2.5) * r
    end
    smap_r /= a

    (isnan(smap_r) || isnan(smap_z)) && throw("nan")

    smap_x = smap_r * (r == 0 ? 0 : x / r)
    smap_y = smap_r * (r == 0 ? 0 : y / r)

#    phi = atan2(y, x)
#    smap_x = smap_r * cos(phi)
#    smap_y = smap_r * sin(phi)

    return Float32(smap_x), Float32(smap_y), Float32(smap_z)
end
