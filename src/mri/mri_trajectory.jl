#=
mri_trajectory.jl
based on mri_trajectory.m
Jing Dong
2020-03
=#

export mri_trajectory

#using MIRT: image_geom_mri, todo
using Random: seed!
using Plots

"""
    kspace, omega, wi = mri_trajectory(ktype, arg_traj, N, fov, arg_wi)

Generate kspace trajectory samples and density compensation functions.
in
`ktype		symbol	k-space trajectory type.  for proper display with help.
`arg_traj	cell	arguments for a specific trajectory
`N	[1 2|3]		target image size
`fov	[1 2|3]		field of view in x and y (and z)
`arg_wi		cell	options to pass to `ir_mri_density_comp`

out
`kspace	[Nk 2|3]	kspace samples in units 1/fov
`omega	[Nk 2|3]	trajectory samples over [-pi,pi]
`wi	[Nk 1]		(optional) density compensation factors
trajectory types:
`:cartesian :radial :cart_y_2 :random
:half_8 :epi_sin
:spiral0 :spiral1 :spiral3
:rosette3
:epi_under
:gads` (emulate golden-angle data sharing per winkelmann:07:aor)
"""
function mri_trajectory(wi ; ktype::Symbol, N::Dims, fov,
  arg_wi, arg_traj, samp::Array = [], na_nr::Real = 2*pi, na::Array = [],
  nr::Real = maximum(N)/2, ir::Array = [],
  omax::Real = pi, Nro::Int = -1,
    delta_ro::Real = 1/Nro,
    shift::Real = -0.75,
    kmax_frac::Array = [0.20 0.35 0.501],
    nspoke::Array = [],
    under::Array = [1 1 0.6])

  if Nro == -1
    Nro = maximum(collect(N))
  end

  if nspoke == []
    nspoke = pi * kmax_frac * Nro
    for i in eachindex(nspoke)
      nspoke[i] = floor(nspoke[i])
    end
  end

    if (length(N) == 1)
        N = [N N]
    end

    if (length(fov) == 1)
        fov = fov * ones(size(N))
    end

    #default density compensation, which works for ordinary Cartesian (only)
    #wi = []

    #=
    trajectory choices
    ideal cartesian
    =#
    if (ktype == :cartesian)

       #todo: antonis will fix
       o1 = (collect(0:(N[1]-1))/N[1] .- 0.5)*2*pi
       o2 = (collect(0:(N[2]-1))/N[2] .- 0.5)*2*pi
       if (length(N) == 2)
         o1, o2 = ndgrid(o1, o2)
         omega = [o1[:] o2[:]]
       end

       if (length(N) == 3)
         o3 = (0:(N[3]-1)/N[3] - 0.5)*2*pi
         o1, o2, o3 = ndgrid(o1, o2, o3)
         omega = [o1[:] o2[:] o3[:]]

       else
         @warn("fail:only 2d and 3d done")
       end

       wi = 1 / prod(fov)
    end

    if (ktype == :epi_under)
      omega, wi = mri_trajectory_epi_under(N, fov, samp = samp)
    end

    if (ktype == :gads)
       omega, wi = mri_trajectory_gads(N, fov, Nro = Nro,
       delta_ro = delta_ro, shift = shift, kmax_frac = kmax_frac,
       nspoke = nspoke, under = under)
     end

    if (ktype == :radial)
       omega, wi = mri_trajectory_radial(N = N, fov = fov, na_nr = na_nr,
       na = na, nr = nr, ir = ir, omax = omax)

       if (length(N) == 3) # 'stack of radials'
         omega = mri_trajectory_stack(omega, N[3])
         wi = repeat(wi, N[3], 1)
         wi = wi[:] / fov[3]
       end

       #wi = [] # discard the default "analytical" DCFs

     end

    if (ktype == :rosette3)
       if length(N) != 3
         @warn("fail:only 3d done")
       end
       omega, wi = mri_trajectory_rosette3(N, fov)
     end

    #half cartesian + 8 rows
    if (ktype == :half_8)
       o1 = (collect(0:(N[1]-1))/N[1] .- 0.5)*2*pi
       o2 = collect(-N[2]/2:8)/N[2] * 2*pi
       omega = ndgrid(o1, o2)
       @show(typeof(omega))
    end

    # echo-planar with sinusoid:
    if (ktype == :epi_sin)
      if (length(arg_traj) == 0)
         oversample = 1

      elseif (length(arg_traj) == 1) #omit the iscell in matlab code
         oversample = arg_traj[1]
      end
       Npt = oversample*prod(N)
       t = collect(0:(Npt-1))/Npt
       omega = [pi*sin.(2*pi*t*N[2]/2) t*2*pi.-pi]
       @show(size(omega))
     end

    # bad spiral:
    if (ktype == :spiral0)
       Nspiral = 400
       omega = transpose(range(0, stop = 10*2*pi, length = Nspiral))
       omega = pi*[cos.(omega) sin.(omega)] .* [omega omega]/maximum(omega)
       if isempty(arg_wi)
         wi = abs.(omega[:,1] + 1im * omega[:,2]) # simple |r| weighting
       end
     end

    # crude spiral:
    if (ktype == :spiral1)
       Nspiral = round(prod(N) * pi/4)
       omega = range(0, stop = N[1]*2*pi, length = convert(UInt128,Nspiral))
       #@show(size(omega))
       max_omega = maximum(omega)
       #@show(max_omega)
       omega = pi*[cos.(omega) sin.(omega)].* [omega omega]/max_omega
       #@show(omega)
       if isempty(arg_wi)
         wi = abs.(omega[:,1] + 1im * omega[:,2]) # simple |r| weighting
       end
    end

    # 3T spiral:
    if (ktype == :spiral3)
       kspace, omega = mri_kspace_spiral(N = Int(maximum(N[1:2])), fov = maximum(fov[1:2]))

       if length(N) == 3 #stack of these spirals
         omega = mri_trajectory_stack(omega, N[3]);
       end

       wi = abs.(omega[:,1] + 1im * omega[:,2]); # simple |r| weighting
       wi = wi / (fov[1]*fov[2]); # approximate scaling
       if length(N) == 3
         wi = wi / fov[3]; # cartesian-style weighting in z
       end
     end

    #random
    if (ktype == :random)
       seed!(0)
       omega = (rand(N[1]*N[2]*2, 2).-0.5)*2*pi;
    end

    #2D FT, undersampled in "y" (phase encode) direction
    if (ktype == :cart_y_2)
       o1 = (collect(0:(N[1]/1-1))/(N[1]/1) .- 0.5)*2*pi
       o2 = (collect(0:(N[2]/2-1))/(N[2]/2) .- 0.5)*2*pi
       omega = ndgrid(o1, o2)
    end

    # convert to physical units
    if (ktype == :half_8 || ktype == :cart_y_2)
      dx = fov[1] / N[1]
      kspace = omega ./ (2*pi) ./ dx
    else
      kspace = zeros(size(omega))
      #for id = 1:length(N)
       id = 1
       dx = fov[id] / N[id]
       kspace[:,id] = omega[:,id] / (2*pi) / dx
    end

    return kspace, omega, wi
end


# mri_trajectory_stack()
# make 3D "stack of ..." trajectory from a 2D trajectory
function mri_trajectory_stack(omega2, N3::Real)
  o3 = (collect(0:(N3-1)/N3) .- 0.5)*2*π
  o3 = repeat(o3, sizeof(omega2)*2, 1) # [N12,N3]
  omega = repeat(omega2, N3, 1) # [N12*N3,2]
  omega = [omega o3] # [N12*N3,3]
return omega
end


# mri_trajectory_epi_under()
# EPI, with optional under-sampling
function mri_trajectory_epi_under(N, fov; samp::Array = [])
  # samp = true(N(2),1) default keeps all phase-encode samples
  if(samp == [])
    samp = trues(N[2],1)
  end
  nx = N[1]
  ny = N[2]
  omx = collect(-nx/2:nx/2-1) / nx * 2*pi # [-pi,pi) in x
  omy = (ny-1-ny/2) / ny * 2 * pi
  @show(length(omx))
  @show(length(omy))
  omega = [omx omy*ones(nx,1)]
  for iy = 1:ny  # of possible phase encodes
   if samp[iy]
     omy = (iy-1-ny/2) / ny * 2 * pi
     omega = [omega; [omx omy*ones(nx,1)]]
     omx = reverse(omx, dims = 1)
   end
  end
  wi = ones(size(omega,1)) / prod(fov)
return omega, wi
end


"""
    mri_trajectory_gads()
emulate 2D golden angle radial sampling with data sharing
Arguments:
Nro: of samples in each readout/spoke
shift: shift along read-out due to gradient delays (stress)
kmax_frac: fractions of maximum krad (0.5) for rings
under: under-sampling factor for each annulus
Output:
omega, wi
"""
function mri_trajectory_gads(N, fov; Nro::Int = -1,
  delta_ro::Real = 1/Nro,
  shift::Real = -0.75,
  kmax_frac::Array = [0.20 0.35 0.501],
  nspoke::Array = floor(pi * kmax_frac * Nro),
  under::Array = [1 1 0.6])
  delta_ro = 1/Nro
  if (Nro == -1)
    Nro = maximum(collect(N))
  end
  nspoke = nspoke .* under
  if true # make fibonacci for more uniform coverage
   #@info(nspoke)
   phi = (1 + sqrt(5)) / 2;
   n = round.(log.(nspoke * sqrt(5) .+ 0.5) / log(phi))
   nspoke = phi.^n / sqrt(5) .+ 0.5
   for i in eachindex(nspoke)
     nspoke[i] = floor(nspoke[i])
   end
  end
  start = [0 1 2] * pi/4
  nring = length(nspoke)
  kmax_frac = [0 kmax_frac]
  omega = zeros(0,2)
  for ir=1:nring
    #TODO: not sure here
     kspace = ir_mri_kspace_ga_radial(Nspoke = Int(nspoke[ir]),
     Nro = Nro, delta_ro = delta_ro, shift = shift,
       start = start[ir])

     kspace = reshape(kspace, :, 2)
     krad = sqrt.(sum(kspace.^2, dims = 2)[:])
     good = (kmax_frac[ir] .<= krad) .& (krad .< kmax_frac[ir+1])
     kspace = kspace[good, :]
     omega = [omega; 2*pi*kspace]
  end
  wi = [] # no default DCF
  return omega, wi
end


"""
    mri_trajectory_radial()
Arguments:
# na_nr: default ensures proper sampling at edge of k-space
# na: angular spokes (default: na_nr * nr)
# nr: radial samples per spoke
# ir: default: 0:nr
# omax: maximum omega
# todo: generalize to 3D using barger:02:trc
"""
function mri_trajectory_radial(;N, fov,
  na_nr::Real = 2*pi, na::Array = [],
  nr::Real = maximum(N)/2, ir::Array = [],
  omax::Real = pi)
  if isempty(ir)
    ir = collect(0:nr)
  end
  if isempty(na)
    na = 4*ceil(na_nr * nr/4) # mult of 4
  end
  om = ir/nr * pi
  ang = collect(0:na-1)/na * 2*pi
  om,ang = ndgrid(om, ang) # [nr+1, na]
  #@show(size(reshape(om.*cos.(ang),1716)))
  omega = [reshape(om.*cos.(ang),1716) reshape(om.*sin.(ang),1716)]

# density compensation factors based on "analytical" voronoi

  du = 1/fov[1] # assume this radial sample spacing
  wi = pi * du^2 / na * 2 * ir[:] # see lauzon:96:eop, joseph:98:sei
  for i in collect(1:length((ir)))
    if ir[i] == 0
      wi[i] = pi * (du/2)^2 / na
    end
  end
  return omega, wi
end

  #wi = repeat(wi, convert(UInt128,na))
  #@show(size(wi))
  #wi = wi[:]


"""
   mri_trajectory_rosette3()
3d rosette, with default parameters from bucholz:08:miw
Arguments:
omax: maximum omega
nt : time samples (65.536 ms for 4 usec dt)
dt : time sample spacing (4 usec)
ti : time samples
"""
function mri_trajectory_rosette3(N, fov; f1::Int = 211,
  f2::Real = 117.13, f3::Real = 73.65,
  nshot::Int = 32, omax::Real = pi, nt::Int = 16385,
  dt::Real = 4e-6, ti::Array = [])
  #arg = vararg_pair(arg, varargin);
  ti = collect(0:nt-1) * dt

  tmp = 2 * pi * ti;
  p1 = f1 * tmp;
  p2 = f2 * tmp;
  p3 = f3 * tmp;
  kx = omax * sin.(p1) .* cos.(p2) .* cos.(p3)
  ky = omax * sin.(p1) .* sin.(p2) .* cos.(p3)
  kz = omax * sin.(p1) .* sin.(p3)
  omega = [kx ky kz]
  for is= collect(1:nshot-1) # n-shot, rotate kx,ky by 2 pi / N
   ang = is * 2 * pi / nshot;
   c = cos(ang);
   s = sin(ang);
   ox = c * kx + s * ky;
   oy = -s * kx + c * ky;
   omega = [omega; [ox oy kz]];
  end
  wi = omax^3 * abs.(sin.(p1).^2 .* cos.(p1) .* cos.(p3)) #from bucholz:08:miw
  return omega, wi
end


"""
   mri_trajectory(:test)
self test
"""
function mri_trajectory(test::Symbol)
  test != :test && throw(DomainError(test, "Not valid"))
  ig = image_geom_mri(nx = 2^6, ny = 2^6-0, fov = 250) # 250 mm FOV
  N = ig.dim

  arg_tr = []
  arg_wi = []
  ktype = [:spiral3,:cartesian, :radial, :gads, :rosette3, :epi_sin,
            :spiral0, :spiral1, :random, :epi_under]
  #arg_tr = {:na_nr, pi/2}
  arg_wi = [:voronoi]
  wi = []
  plots = Any[]
  for i in 1:length(ktype)
    kspace, omega, wi = mri_trajectory(arg_tr, ktype = ktype[i],
       N = N, fov = ig.fovs, arg_wi = arg_wi, arg_traj = arg_tr, na_nr = pi/2)
    @show(ktype[i])
    push!(plots, plot(omega[:,1], omega[:,2],
            xlabel = "omega1",
            ylabel = "omega2", label = string(ktype[i])))

    if (ktype[i] == :epi_sin)
      #plot wi != 0case
      arg_tr = [10]
      kspace, omega, wi = mri_trajectory(arg_tr, ktype = ktype[i],
      N = N, fov = ig.fovs, arg_wi = arg_wi, arg_traj = arg_tr, na_nr = pi/2)
      push!(plots, plot(omega[:,1], omega[:,2],
              xlabel = "omega1",
              ylabel = "omega2", label = string(ktype[i])))
    end
  end
  plot(plots...)

  ktype = [:half_8, :cart_y_2]
  plots = Any[]
  for i in 1:length(ktype)
    kspace, omega, wi = mri_trajectory(arg_tr, ktype = ktype[i],
    N = N, fov = ig.fovs, arg_wi = arg_wi, arg_traj = arg_tr, na_nr = pi/2)
    @show(ktype[i])
    omega1 = omega[1]
    omega2 = omega[2]
    push!(plots, plot(omega1[:], omega2[:],
          xlabel = "omega1",
          ylabel = "omega2", label = string(ktype[i])))
  end
  plot(plots...)

  arg_wi = []
  ktype = [:spiral0, :spiral1]
  plots = Any[]
  for i in 1:length(ktype)
    kspace, omega, wi = mri_trajectory(arg_tr, ktype = ktype[i],
    N = N, fov = ig.fovs, arg_wi = arg_wi, arg_traj = arg_tr, na_nr = pi/2)
    @show(ktype[i])
    push!(plots, plot(omega[:,1], omega[:,2],
            xlabel = "omega1",
            ylabel = "omega2", label = string(ktype[i])))
  end
  plot(plots...)

  ig = image_geom_mri(nz = 2^5, nx = 2^6, ny = 2^6-0, fov = 250, dz = 2^5) # 250 mm FOV
  N = ig.dim
  ktype = [:cartesian, :radial, :gads, :spiral3]
  plots = Any[]
  for i in 1:length(ktype)
    kspace, omega, wi = mri_trajectory(arg_tr, ktype = ktype[i],
    N = N, fov = ig.fovs, arg_wi = arg_wi, arg_traj = arg_tr, na_nr = pi/2)
    @show(ktype[i])
    push!(plots, plot(omega[:,1], omega[:,2],
            xlabel = "omega1",
            ylabel = "omega2", label = string(ktype[i])))
  end
  plot(plots...)
  return true
end
