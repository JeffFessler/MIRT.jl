#=
kspace_spiral.jl
Jing Dong

Based on mri_kspace_spiral.m that was based on m-files from Valur Olafsson
that he got from Brad Sutton who got them from Doug Noll...
=#

export mri_kspace_spiral

using Interpolations


"""
    kspace, omega, gradxy = mri_kspace_spiral( [options] )

Make k-space spiral trajectory based on GE 3T scanner constraints

Option:
- `N` dimention of reconstructed image
- `Nt` # of time points
- `fov` field of view in cm
- `dt` time sampling interval out; default `5e-6` sec
- `gamp::Real` design gradient amplitude max, G/cm; default 2.2
- `gslew::Int` design slew rate, mT/m/ms; default 180

Out:
- `kspace [Nt,2]` kspace trajectory `[kx ky]` in cycles/cm, NO: cycles/FOV
- `omega [Nt,2]` "" in radians
- `gradxy [Nt 2]` gradient waveforms in (units?)
"""
function mri_kspace_spiral( ;
    fov::Real = 22,
    N::Int = 64, Nt::Int = -1,
    dt::Float64 = 5e-6, nl::Int = 1,
    gamp::Real = 2.2, gslew::Int = 180,
    warn_nk::Bool = true,
)

    if Nt == -1
        if fov == 20
            Nt = 4026
        elseif fov == 22
            Nt = 3770
        else
            Nt = 0 # let algorithm choose
        end
    end

    # generate spiral k-space trajectory
    kx, ky, gx, gy = genkspace(fov, N, Nt, nl, gamp, gslew, dt, warn_nk)
    if (nl == 1)
        kspace = [kx ky]
        gradxy = [gx gy]
    else
        kspace = permutedims(cat(dims = 3, kx, ky), [1,3,2])
        gradxy = permutedims(cat(dims = 3, gx, gy), [1,3,2])
    end
    omega = 2π * [kx ky] / N

    maximum(omega) > π && throw("bad spiral")

    return kspace, omega, gradxy
end


"""
    genkspace

Generate the proper length of k-space trajectory.

It linearly interpolates the output of `genspiral` to the correct `length()`
& takes care of the rotations for the interleaves.
* `ld` is the length of the data
* `nint` is the number of interleaves
Brad Sutton; University of Michigan
"""
function genkspace(
    FOV, N, ld, nint, gamp, gslew, tsamp, warn_nk ;
    rotamount::Int = 0,
)
    nk = ld/nint

    flag = (nk == 0) # auto determine number of k-space points
    # just input ld = 0.3

    #dts = 4e-6 # 5e-6 [sec]

    Gx, Gy, kxi, kyi, sx, sy, dts =
        genspi(FOV, N, nl=nint, gamp=gamp, gslew=gslew)

    ik = 0:(length(kxi)-1)
    tk = 0:(dts/tsamp*length(kxi)-1)
    tk = tk * tsamp
    kxt = interp1(ik*dts, kxi, tk)
    kyt = interp1(ik*dts, kyi, tk)

    ig = 0:(length(Gx)-1)
    tg = 0:(dts/tsamp*length(Gx)-1)
    tg = tg * tsamp
    gxt = interp1(ig*dts, Gx, tg)
    gyt = interp1(ig*dts, Gy, tg)

    if flag == 1
        nk = length(kxt) - 2
    end
    nk = Int64.(nk)

    if (nk > length(kxt))
        warn_nk && @warn "reduce nk from $nk to $(length(kxt))"
        nk = min(nk, length(kxt))
    end
    kxo = kxt[1:nk]
    kyo = kyt[1:nk]

    gxo = gxt[1:nk]
    gyo = gyt[1:nk]

    # rotate matrix for proper orientation
    phir = -rotamount*pi/2
    kxop = kxt*cos(phir) - kyt*sin(phir)
    kyop = kyt*cos(phir) + kxt*sin(phir)
    gxop = gxt*cos(phir) - gyt*sin(phir)
    gyop = gyt*cos(phir) + gxt*sin(phir)

    kx = zeros(nk, nint)
    ky = zeros(nk, nint)
    gx = zeros(nk, nint)
    gy = zeros(nk, nint)

    if (length(kxop) > length(nk))
        kx = zeros(length(kxop), nint)
        ky = zeros(length(kyop), nint)
        gx = zeros(length(gxop), nint)
        gy = zeros(length(gyop), nint)
    end

    kx[:,1] = kxop
    ky[:,1] = kyop
    gx[:,1] = gxop
    gy[:,1] = gyop
    phi = 2*π/nint
    if nint > 1
        phi = 2*π/nint
        for ii = 1:(nint-1)
            kx[:,ii+1] = kxop*cos(ii*phi) - kyop*sin(ii*phi)
            ky[:,ii+1] = kyop*cos(ii*phi) + kxop*sin(ii*phi)
            gx[:,ii+1] = gxop*cos(ii*phi) - gyop*sin(ii*phi)
            gy[:,ii+1] = gyop*cos(ii*phi) + gxop*sin(ii*phi)
        end
    end

    return kx, ky, gx, gy
end


"""
    Gx, Gy, kx, ky, sx, sy, gts = genspi(...)

This is translation of C code from scanner:
exactly what is played out
to gradients at 4us.

multi-shot spiral design
uses Duyn's approximate slewrate limited design
augmented with archimedian `gmax` limit
in [args]
* `D` = FOV; cm
* `N` = matrix size()

* `Tmax` = longest acquisition allowed; s
* `dts` = output sample spacing; s
* `gtype` = trajectory type()

option [CVs]
* `nl` = number of interleaves
* `gamp` = design grad max; G/cm
* `gslew` = design slew rate; mT/m/ms
* `nramp` = number of rampdown points; default 0

out
* `Gx; Gy`

time is in sec()
* rev 0 12/26/98    original
* rev 1 4/15/99    little better calc of ts

Borrowed from Doug Noll; Univ. of Michigan.
Modified to take more input cv's.
"""
function genspi(D, N ;
    nl::Int = 1, gamp::Real = 202, gslew::Int = 180, nramp::Int = 0)

    ########## Predefined variables
    GRESMAX = 21000
    # nramp=100
    MAX_PG_WAMP = 32766
    gts = 4e-6 # [sec]

    Tmax = GRESMAX * gts
    dts = gts
    opfov = D

    #################################
    gamma = 2*π*4.257e3
    gambar = gamma / (2*π)

    gx = zeros(2*GRESMAX)
    gy = zeros(2*GRESMAX)

    q = 5
    S0 = gslew*100
    dt = dts*.5

    # slew-rate limited approximation

    Ts = 0.666667 / nl*sqrt(((pi*N)^3)/(gamma*D*S0))

    a2 = N*π/(nl*(Ts^(0.666667)))
    a1 = 1.5*S0/a2
    beta = S0*gamma*D/nl
    Gmax = a1*(Ts^0.333333)
    gmax = 0

    t = 0:dt:Ts
    x = t .^ 1.333333
    theta = (t.^2) .* (.5 * beta ./ (q .+ .5 * beta ./ a2 .* x))
    y = q .+ .5 .* beta./a2 .* x
    dthdt = t .* (beta .* (q .+ 0.166667 * beta ./ a2 .* x) ./ (y .* y))
    c = cos.(theta)
    s = sin.(theta)
    gx = (nl/(D * gamma)) .* dthdt .* (c - theta .* s)
    gy = (nl/(D * gamma)) .* dthdt .* (s + theta .* c)
    gabs = abs.(gx + 1im*gy)
    # cut short if over peak
    gmax = abs.(gamp ./ (theta .+ eps()).+ 1im * gamp)
    l1 = length(t) - sum(gabs .> gmax)
    ts = t[l1]
    thetas = theta[l1]

    # gmax limited approximation
    l3 = 0
    T = ts
    if Gmax > gamp
        T = ((π*N/nl)*(π*N/nl) - thetas*thetas)/(2*gamma*gamp*D/nl) + ts
        t = (ts+dt):dt:T
        theta = sqrt.(thetas * thetas .+ (2 * gamma * gamp * D) .* (t .- ts) / nl)
        c = cos.(theta)
        s = sin.(theta)
        ind2 = l1 .+ (1:length(t))
        if (maximum(ind2) > length(gx))
            gx = [gx; zeros(maximum(ind2) - length(gx))]
        end
        if (maximum(ind2) > length(gy))
            gy = [gy; zeros(maximum(ind2) - length(gy))]
        end
#       gx[floor.(Int,ind2)] = gamp.*(c./theta - s)
#       gy[floor.(Int,ind2)] = gamp.*(s./theta + c)
        gx[ind2] = gamp.*(c./theta - s)
        gy[ind2] = gamp.*(s./theta + c)
        l3 = length(t)
    end

    l2 = l1 + l3
    Gx = gx[1:2:l2] # | gx[1:2:l2]*MAX_PG_WAMP/gamp
    Gy = gy[1:2:l2] # | gy[1:2:l2]*MAX_PG_WAMP/gamp
    g = Gx + 1im *Gy # slew rate vector
    s = diff(g)./(gts*1000)  # grad vector
    Kx = cumsum(Gx)*gts*opfov*gambar
    Ky = cumsum(Gy)*gts*opfov*gambar
    k = Kx + 1im*Ky  # kspace vector
    t = collect(0:gts:T) # time vector
    matrix = maximum(abs.(k))*2
    maxg = maximum(abs.(g))
    maxs = maximum(abs.(s))
    maxt = maximum(t)*1000

    kx = real(k)
    ky = imag(k)
    sx = real(s)
    sy = imag(s)

    return Gx, Gy, kx, ky, sx, sy, gts
end
