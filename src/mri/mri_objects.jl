#=
mri_objects.jl
2019-11 Connor Martin
based on mri_objects.m
Copyright 2007-6-28, Jeff Fessler, University of Michigan
=#

export mri_objects

using Plots
using FFTW: fft, fftshift, ifftshift
using MIRT: max_percent_diff,jinc,rect
#using MIRT: jim, image_geom, prompt, interp1
#include("max_percent_diff.jl")
#include("jinc.jl") # todo - put in MIRT
#include("rect.jl")

RealArray = AbstractArray{<:Real}
PairList = Vector{<:Tuple} # Vector{Tuple{Symbol,<:Any}}

# struct for wrapper definition
struct MIRT_mri_object
    pairs::PairList
    image::Function
    kspace::Function
end


# the two primary functions
# nvp is a vector of tuples (type::Symbol, parameter vector)
# for 2d, just make the third coord 0

function mri_objects_image(nvp::PairList, x::RealArray, y::RealArray, z ;
        dx::Real=0, dy::Real=0, dz::Real=0)
    out = zeros(size(x)...)
    for pair in nvp
        out += imageFuns[pair[1]](pair[2], x, y, z ; dx=dx, dy=dy, dz=dz)
    end
    return out
end

function mri_objects_kspace(nvp::PairList, u::RealArray, v::RealArray, w)
    out = zeros(size(u)...)
    for pair in nvp
        out += kspaceFuns[pair[1]](pair[2], u, v, w)
    end
    return out
end

# 2D
function mri_objects_image(data, x::RealArray, y::RealArray ; kwargs...)
    return mri_objects_image(data, x, y, 0 ; kwargs...)
end

function mri_objects_kspace(data, u::RealArray, v::RealArray)
    return mri_objects_kspace(data, u, v, 0)
end


# The big chunk of functions for each type (image & kspace)

function mri_objects_trap(z, dz, len)
    dz = abs(dz)
    (dz == 0 || len < dz) && return rect.(z / len)
    interp1([-len/2-dz/2, -len/2+dz/2, len/2-dz/2, len/2+dz/2], [0, 1, 1, 0], z)
end


function mri_objects_image_gauss2(data::RealArray, x::RealArray, y::RealArray,
        z ; kwargs...)
    size(data,2) != 5 && throw(DomainError(data,"gauss2 requires 5 parameters"))
    z = zeros(size(data,1))
    zw = z .+ Inf
    inputParams = [data[:,1:2] z data[:,3:4] zw data[:,5]]
    return mri_objects_image_gauss3(inputParams,x,y,0)
end

function mri_objects_image_gauss3(data::RealArray, x::RealArray, y::RealArray,
        z ; kwargs...)
    size(data,2) != 7 && throw(DomainError(data,"gauss3 requires 7 parameters"))
    out = zeros(size(x))
    for i=1:size(data,1)
        par = data[i,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        xw = par[4] ./ sqrt(log(256)) .* sqrt(2*pi)
        yw = par[5] ./ sqrt(log(256)) .* sqrt(2*pi)
        zw = par[6] ./ sqrt(log(256)) .* sqrt(2*pi)
        param1 = exp.(-pi .* ((x .- xc) ./ xw).^2)
        param2 = exp.(-pi .* ((y .- yc) ./ yw).^2)
        param3 = exp.(-pi .* ((z .- zc) ./ zw).^2)
        out = out .+ par[7] .* param1 .* param2 .* param3
    end
    return out
end

function mri_objects_image_circ2(data::RealArray, x::RealArray, y::RealArray,
        zr ; kwargs...)
    size(data,2) != 4 && throw(DomainError(data,"circ2 requires 4 parameters"))
    z = zeros(size(data,1))
    inputParams = [data[:,1:2] z data[:,3] (z .+ 1) data[:,4]]
    return mri_objects_image_cyl3(inputParams, x,y,zr ; kwargs...)
end

function mri_objects_image_cyl3(data::RealArray, x::RealArray, y::RealArray,
        z ; dx::Real=0, dy::Real=0, dz::Real=0)
    size(data,2) != 6 && throw(DomainError(data,"cyl3 requires 6 parameters"))
    out = zeros(size(x))
    for i = 1:size(data,1)
        par = data[i,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        rad = par[4]
        len = par[5]
        out += par[6] * ((x .- xc).^2 + (y .- yc).^2 .< rad^2) .*
            mri_objects_trap(z .- zc, dz, len)
    end
    return out
end

function mri_objects_image_rect2(data::RealArray, x::RealArray, y::RealArray,
        zr ; dx::Real=0, dy::Real=0, dz::Real=0)
    size(data,2) != 5 && throw(DomainError(data,"rect2 requires 5 parameters"))
    z = zeros(size(data,1))
    inputParams = [data[:,1:2] z data[:,3:4] (z .+ 1) data[:,5]]
    return mri_objects_image_rect3(inputParams, x,y,zr ; dx=dx, dy=dy, dz=dz)
end

function mri_objects_image_rect3(data::RealArray, x::RealArray, y::RealArray,
        z ; dx::Real=0, dy::Real=0, dz::Real=0)
    out = zeros(size(x))
    for i=1:size(data,1)
        par = data[i,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        xw = par[4]
        yw = par[5]
        zw = par[6]
        out = out .+ par[7] *
            mri_objects_trap(x .- xc, dx, xw) .*
            mri_objects_trap(y .- yc, dy, yw) .*
            mri_objects_trap(z .- zc, dz, zw)
    end
    return out
end

function mri_objects_image_dirac2(data::RealArray, x::RealArray, y::RealArray,
        z ; kwargs...)
    size(data,2) != 3 && throw(DomainError(data,"dirac2 requires 3 parameters"))
    #@show data
    inputParams = [data[:,1:2] zeros(size(data,1),1) data[:,3]]
    return mri_objects_image_dirac3(inputParams, x,y,z ; kwargs...)
end

function mri_objects_image_dirac3(data::RealArray, x::RealArray, y::RealArray,
        z ; kwargs...)
    size(data,2) != 4 && throw(DomainError(data,"dirac3 requires 4 parameters"))
    out = zeros(size(x))
    for i=1:size(data,1)
        out += data[i,4] .* ((x .== data[i,1]) .* (y .== data[i,2]) .* (z .== data[i,3]))
    end
#TODO: what does out(out != 0) = inf, warn do? Ask fess. todo: test
    replace!(x -> x==0 ? 0 : Inf, out) #no temp array
    return out
end


function mri_objects_kspace_gauss2(data::RealArray, u::RealArray, v::RealArray, w)
    size(data,2) != 5 && throw(DomainError(data,"gauss2 requires 5 parameters"))
    z = zeros(size(data,1))
    zw = (1 .+ z) .* sqrt(log(256)) ./ sqrt(2*pi)
    inputParams = [data[:,1:2] z data[:,3:4] zw data[:,5]]
    return mri_objects_kspace_gauss3(inputParams,u,v,w)
end

function mri_objects_kspace_gauss3(data::RealArray, u::RealArray, v::RealArray, w)
    size(data,2) != 7 && throw(DomainError(data,"gauss3 requires 7 parameters"))
    out = 0
    for ii=1:size(data,1)
        par = data[ii,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        xw = par[4] / sqrt(log(256)) * sqrt(2*pi)
        yw = par[5] / sqrt(log(256)) * sqrt(2*pi)
        zw = par[6] / sqrt(log(256)) * sqrt(2*pi)
        out = out .+ par[7] * xw * yw * zw *
            exp.(-pi * (u .* xw) .^2) .*
            exp.(-pi * (v .* yw) .^2) .*
            exp.(-pi * (w .* zw) .^2) .*
            exp.(-2im * pi * (u .* xc .+ v .* yc .+ w .* zc))
    end
    return out
end

function mri_objects_kspace_circ2(data::RealArray, u::RealArray, v::RealArray, w)
    size(data,2) != 4 && throw(DomainError(data,"circ2 requires 4 parameters"))
    z = zeros(size(data,1))
    inputParams = [data[:,1:2] z data[:,3] (z .+ 1) data[:,4]]
    return mri_objects_kspace_cyl3(inputParams,u,v,w)
end

function mri_objects_kspace_cyl3(data::RealArray, u::RealArray, v::RealArray, w)
    size(data,2) != 6 && throw(DomainError(data,"cyl3 requires 6 parameters"))
    out = zeros(size(u))
    for i = 1:size(data,1)
        par = data[i,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        rad = par[4]
        len = par[5]
        out = out .+ par[6] * rad^2 * 4 .*
            jinc.(2 * sqrt.(u.^2 .+ v.^2) .* rad) .* (len .* sinc.(w .* len)) .*
            exp.(-2im * pi * (u .* xc .+ v .* yc .+ w .* zc))
    end
    return out
end

function mri_objects_kspace_rect2(data::RealArray, u::RealArray, v::RealArray, w)
    size(data,2) != 5 && throw(DomainError(data,"rect2 requires 5 parameters"))
    z = zeros(size(data,1),1)
    inputParams = [data[:,1:2] z data[:,3:4] (z .+ 1) data[:,5]]
    return mri_objects_kspace_rect3(inputParams,u,v,0)
end

function mri_objects_kspace_rect3(data::RealArray, u::RealArray, v::RealArray, w)
    size(data,2) != 7 && throw(DomainError(data,"rect3 requires 7 parameters"))
    out = 0
    for ii=1:size(data,1)
        par = data[ii,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        xw = par[4]
        yw = par[5]
        zw = par[6]
        out = out .+ par[7] * xw * yw * zw *
            sinc.(u .* xw) .* sinc.(v .* yw) .* sinc.(w .* zw) .*
            exp.(-2im * pi * (u .* xc + v .* yc .+ w .* zc))
    end
    return out
end

function mri_objects_kspace_dirac2(data::RealArray, u::RealArray, v::RealArray, w)
    size(data,2) != 3 && throw(DomainError(data,"dirac2 requires 3 parameters"))
    inputParams = [data[:,1:2] zeros(size(data,1),1) data[:,3]]
    return mri_objects_kspace_dirac3(inputParams,u,v,w)
end

function mri_objects_kspace_dirac3(data::RealArray, u::RealArray, v::RealArray, w)
    size(data,2) != 4 && throw(DomainError(data,"dirac3 requires 4 parameters"))
    out = zeros(size(u))
    for i=1:size(data,1)
        out += data[i,4] * exp.(-2im * pi * (u .* data[i,1] .+ v .* data[i,2] .+ w .* data[i,3]))
    end
    return out
end


function mri_objects_case1( ; unit::Symbol = :mm)
    rp = [ # rect2
     0 0    200 200  1;
    -50 -50  40  40  1;
     50 -50  20  20  1;
     0 50    50  50  1;
    ]
    gp = [ # gauss2
     -70 0  1 1  1;
     -60 0  2 2  1;
     -50 0  3 3  1;
     -40 0  4 4  1;
     -20 0  5 5  1;
       0 0  6 6  1;
      20 0  7 7  1;
      50 0  8 8  1;
    ]

    if unit === :cm
        rp[:,1:6] ./= 10
        gp[:,1:6] ./= 10
    end
    return [(:rect2,rp), (:gauss2,gp)] # PairList
end


# image function dictionary
imageFuns = Dict([
    (:dirac2, mri_objects_image_dirac2),
    (:dirac3, mri_objects_image_dirac3),
    (:rect2, mri_objects_image_rect2),
    (:rect3, mri_objects_image_rect3),
    (:gauss2, mri_objects_image_gauss2),
    (:gauss3, mri_objects_image_gauss3),
    (:circ2, mri_objects_image_circ2),
    (:cyl3, mri_objects_image_cyl3),
])

# kspace function dictionary
kspaceFuns = Dict([
    (:dirac2, mri_objects_kspace_dirac2),
    (:dirac3, mri_objects_kspace_dirac3),
    (:rect2, mri_objects_kspace_rect2),
    (:rect3, mri_objects_kspace_rect3),
    (:gauss2, mri_objects_kspace_gauss2),
    (:gauss3, mri_objects_kspace_gauss3),
    (:circ2, mri_objects_kspace_circ2),
    (:cyl3, mri_objects_kspace_cyl3),
])



"""
    st = mri_objects([(:type1, params1), (:type2, params2), ...])
Generate object that describes image-domain objects and Fourier domain spectra
of simple structures such as rectangles, disks, superpositions thereof.
These functions are useful for simple "idealized" MRI simulations
where the data is modeled as analytical Fourier samples,
i.e., no field inhomogeneity and no relaxation effects.

in
- `[(type, params), ...]` e.g. `[(:rect2, params), (:gauss3, params), ...]`

type, params:
- `:dirac2 [N 3] [xcent ycent value]`
- `:dirac3 [N 4] [xcent ycent zcent value]`
- `:rect2 [N 5] [xcent ycent xwidth ywidth value]`
- `:rect3 [N 7] [xcent ycent zcent xwidth ywidth zwidth value]`
- `:gauss2 [N 5] [xcent ycent xwidth ywidth value]`
- `:gauss3 [N 7] [xcent ycent zcent xwidth ywidth zwidth value]`
- `:circ2 [N 5] [xcent ycent rad value]`
- `:cyl3 [N 6] [xcent ycent zcent xrad zwidth value]`

All types must be 2D or 3D, not mixed.

out
- `st` struct
* `st.image(x,y)`  returns 2D image-domain (sampled) picture
* `st.image(x,y,z)` same but 3D
* `st.kspace(u,v)` returns 2D Fourier-space samples
* `st.kspace(u,v,w)` same but 3D
"""
function mri_objects(nvp::PairList)
    type2 = (:dirac2, :rect2, :gauss2, :circ2)
    type3 = (:dirac3, :rect3, :gauss3, :cyl3)
    is3 = nvp[1][1] ∈ type3
    if is3
        for pair in nvp
            (pair[1] ∉ type3) && throw("mix of 2d and 3d types")
        end

        return MIRT_mri_object(nvp,
              (x,y,z ; kwargs...) -> mri_objects_image(nvp, x,y,z ; kwargs...),
              (u,v,w ; kwargs...) -> mri_objects_kspace(nvp, u,v,w ; kwargs...))
    else
        for pair in nvp
            (pair[1] ∉ type2) && throw("mix of 2d and 3d types")
        end
        return MIRT_mri_object(nvp,
            (x,y ; kwargs...) -> mri_objects_image(nvp, x,y ; kwargs...),
            (u,v ; kwargs...) -> mri_objects_kspace(nvp, u,v ; kwargs...))
    end
end


"""
    mri_objects((type,params)) for a single tuple
"""
mri_objects(pair::Tuple) = mri_objects([pair,])


# special cases

function mri_objects_case4(fov::AbstractVector{<:Real} ; unit::Symbol=:mm)
    length(fov) != 3 && throw("fov must have length 3")
    cp = [0 0 0 0.4*fov[1] fov[3] 1] # cyl3
    rp = Array{Float16}([ # rect3
      -50 -50  0  40 40 40  1;
       50 -50 40  20 20 50  1;
       0  50 -40  30 30 60  1;
    ])
    rp[:,1:3] = rp[:,1:3]/256 * fov[1]
    rp[:,4:6] = rp[:,4:6]/256 * fov[1]
    gp = Array{Float32}([ # gauss3
       -70 0 0  1 1 1  1;
       -60 0 0  2 2 2  1;
       -50 0 0  3 3 3  1;
       -40 0 0  4 4 4  1;
       -20 0 0  5 5 5  1;
        00 0 0  6 6 6  1;
        20 0 0  7 7 7  1;
        50 0 0  8 8 8  1;
    ])
    gp[:,1:3] = gp[:,1:3]/256 * fov[1]
    gp[:,4:6] = gp[:,4:6]/256 * fov[1]
    if (unit === :cm)
        rp[:,1:6] /= 10
        gp[:,1:6] /= 10
    end
    return [(:cyl3, cp), (:rect3, rp), (:gauss3, gp)] #PairList
end


# 2D test
#NOTE: Dirac will not plot very well. Just make sure it prints a gray image, and that's fine. FFT will also be screwed up.
function mri_objects_test2()
    ig = image_geom(nx = 2^7, ny = 2^7 + 2, dx = 4, offsets = :dsp)

    shift = [0.1, 0.2] .* ig.fovs
    sizes = [.15, .1] .* ig.fovs
    tests = ( # 3 separate tests:
                (:circ2, [shift' sizes[1] 2]),
                (:gauss2, [shift' sizes' 2]),
                (:rect2, [shift' sizes' 2]),
                (:dirac2,[shift' 2])
            )
    ntest = length(tests)
    pl = Matrix{Plots.Plot}(undef, ntest,3)

    for (i,test) = enumerate(tests)
        st = mri_objects(test)
        i2 = st.image(ig.xg, ig.yg, dx=ig.dx, dy=ig.dy)
        s2 = ifftshift(fft(fftshift(i2))) * abs(ig.dx * ig.dy)
        f2 = st.kspace(ig.fg...)
        @show max_percent_diff(f2,s2)
    #   @show round.(maximum.([i2, abs.(f2), abs.(s2)]), digits=2)
        pl[i,1] = jim(i2, title = "$(test[1])")
        pl[i,2] = jim(abs.(s2), title = "|fft|")
        pl[i,3] = jim(abs.(f2), title = "|kspace|")
    end
    plot(pl...)
    return plot(pl...)
end


function mri_objects_test_case1()
    ig = image_geom(nx = 2^7, ny = 2^7 + 2, dx = 4, offsets = :dsp)
    xt = mri_objects(:case1).image(ig.xg, ig.yg)
    jim(xt, title = "case1")
end

function mri_objects_test_case4()
    ig = image_geom(nx = 2^5, ny = 2^5 + 2,nz = 2^5, dx = 1,dz = .5, offsets = :dsp)
    xt = mri_objects(:case4).image(ig.xg, ig.yg, ig.zg)
    jim(xt, title = "case4")
end


# 3D test
#NOTE: Dirac will not plot very well. Just make sure it prints a gray image, and that's fine. FFT will also be screwed up.
function mri_objects_test3()
    ig = image_geom(nx = 2^7, nz = 2^5, offsets = :dsp, dx = 4, dz = 3)
    shift = [0.1, 0.2, 0.3] .* ig.fovs
    sizes = [0.15, 0.1, 0.2] .* ig.fovs
    tests = ((:cyl3, [shift' sizes[1] sizes[3] 2]),
             (:gauss3, [shift' sizes' 2]),
             (:rect3, [shift' sizes' 2]),
             (:dirac3,[shift' 2]))
    ntest = length(tests)
    pl = Matrix{Plots.Plot}(undef, ntest,3)
    for (i,test) in enumerate(tests)
        st = mri_objects(test)
        i3 = st.image(ig.xg,ig.yg,ig.zg, dx=ig.dx, dy=ig.dy, dz=ig.dz)
        s3 = ifftshift(fft(fftshift(i3)))
        s3 *= abs(ig.dx * ig.dy * ig.dz)
        fg = ig.fg
        f3 = st.kspace(fg...)
        @show max_percent_diff(f3,s3)
        pl[i,1] = jim(i3, title = "$(test[1])")
        pl[i,2] = jim(abs.(s3), title = "fft")
        pl[i,3] = jim(abs.(f3), title = "kspace")
        plot(pl[i,:]...)
        prompt()
    end
    return plot(pl...)
end


"""
    mri_objects(key::Symbol ; fov::Real=22, unit::Symbol=:mm)
Special cases

In:
- `fov` Store the following parameter in fov.
- `unit` specify units; default :mm
     Currently, only `:cm` (centimeters) supported as an alternate unit

key choices:
- `:test` run test suite (no other params)
- `:rect2half` centered rectangle of width fov/2
- `:rect3half` 3D version
- `:case1` predefined 2D test case
- `:case4` predefined 3D test case
"""
function mri_objects(key::Symbol ; fov::Real=22, unit::Symbol=:mm)
    if key == :test
        return mri_objects_test()
    elseif key == :case1
        return mri_objects(mri_objects_case1( ; unit=unit))
    elseif key == :case4
        return mri_objects(mri_objects_case4([fov,fov,fov] ; unit=unit))
    elseif key == :rect2half
        return mri_objects((:rect2, [0 0 fov/2 fov/2 1]))
    elseif key == :rect3half
        return mri_objects((:rect3, [0 0 0 fov/2 fov/2 fov/2 1]))
    end
    throw("bad key $key")
end


# plot the trapezoid
function mri_objects_trap_test()
    len = 6
    dz = 2
    z = LinRange(-8,8,201)
    trap = mri_objects_trap(z, dz, len)
    trap0 = mri_objects_trap(z, 0, len)
    plot(z, trap, label="trap")
    plot!(z, trap0, label="rect")
    plot!(xtick = [0 -len/2 len/2-dz/2 len/2 len/2+dz/2])
end


function mri_objects_test()
    @test mri_objects_trap_test() isa Plots.Plot
    prompt()

    @test mri_objects_test2() isa Plots.Plot
    prompt()

	@test mri_objects_test_case1() isa Plots.Plot
    prompt()

    @test mri_objects_test3() isa Plots.Plot

    @test mri_objects_test_case4() isa Plots.Plot
    prompt()
    true
end
