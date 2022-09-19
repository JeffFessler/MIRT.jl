# mri_objects.jl

using MIRT: mri_objects

#using MIRTjim: jim, prompt
using MIRT: max_percent_diff
using ImageGeoms: ImageGeom, fovs, grids, gridf
import MIRT: mri_objects_trap
using FFTW: fft, fftshift, ifftshift
#using LazyGrids: ndgrid
#using Plots
using Test: @test, @test_throws #, @inferred


# 2D test
#NOTE: Dirac will not plot very well.
#Just make sure it prints a gray image, and that's fine. FFT will also be screwed up.
function mri_objects_test2()
    ig = ImageGeom((2^7,2^7 + 2), (4, 4), :dsp)

    shift = [0.1, 0.2] .* fovs(ig)
    sizes = [.15, .1] .* fovs(ig)
    tests = ( # 3 separate tests:
                (:circ2, [shift' sizes[1] 2]),
                (:gauss2, [shift' sizes' 2]),
                (:rect2, [shift' sizes' 2]),
                (:dirac2,[shift' 2])
            )
    ntest = length(tests)
#   pl = Matrix{Plots.Plot}(undef, ntest,3)

    for (i,test) = enumerate(tests)
        st = mri_objects(test)
        i2 = st.image(grids(ig)..., dx=ig.deltas[1], dy=ig.deltas[2])
        s2 = ifftshift(fft(fftshift(i2))) * abs(prod(ig.deltas))
        f2 = st.kspace(gridf(ig)...)
        isinteractive() && (@show max_percent_diff(f2,s2))
    #   @show round.(maximum.([i2, abs.(f2), abs.(s2)]), digits=2)
    #   pl[i,1] = jim(i2, title = "$(test[1])")
    #   pl[i,2] = jim(abs.(s2), title = "|fft|")
    #   pl[i,3] = jim(abs.(f2), title = "|kspace|")
    end
#   plot(pl...)
#   return plot(pl...)
    true
end


function mri_objects_test_case1()
    ig = ImageGeom((2^7,2^7 + 2), (4, 4), :dsp)
    xt = mri_objects(:case1,unit = :cm).image(grids(ig)...)
#   jim(xt, title = "case1cm")
    true
end


function mri_objects_test_case4()
    ig = ImageGeom((2^5,2^5 + 2,2^5), (1, 1, 0.5), :dsp)
    xt = mri_objects(:case4,unit = :cm).image(grids(ig)...)
#   jim(xt, title = "case4cm")
    true
end


# 3D test
# NOTE: Dirac will not plot very well.
# Just make sure it prints a gray image, and that's fine. FFT will also be screwed up.
function mri_objects_test3()
    ig = ImageGeom((2^7,2^7,2^5), (4, 4, 3), :dsp)
    shift = [0.1, 0.2, 0.3] .* fovs(ig)
    sizes = [0.15, 0.1, 0.2] .* fovs(ig)
    tests = ((:cyl3, [shift' sizes[1] sizes[3] 2]),
             (:gauss3, [shift' sizes' 2]),
             (:rect3, [shift' sizes' 2]),
             (:dirac3,[shift' 2]))
    ntest = length(tests)
#   pl = Matrix{Plots.Plot}(undef, ntest,3)
    for (i,test) in enumerate(tests)
        st = mri_objects(test)
        i3 = st.image(grids(ig)..., dx=ig.deltas[1], dy=ig.deltas[2], dz=ig.deltas[3])
        s3 = ifftshift(fft(fftshift(i3))) * abs(prod(ig.deltas))
        f3 = st.kspace(gridf(ig)...)
        isinteractive() && (@show max_percent_diff(f3,s3))
#       pl[i,1] = jim(i3, title = "$(test[1])")
#       pl[i,2] = jim(abs.(s3), title = "fft")
#       pl[i,3] = jim(abs.(f3), title = "kspace")
#       plot(pl[i,:]...)
#       prompt()
    end
#   return plot(pl...)
    true
end


# plot the trapezoid
function mri_objects_trap_test()
    len = 6
    dz = 2
    z = LinRange(-8,8,201)
    trap = mri_objects_trap(z, dz, len)
    trap0 = mri_objects_trap(z, 0, len)
#   plot(z, trap, label="trap")
#   plot!(z, trap0, label="rect")
#   plot!(xtick = [0, -len/2, len/2-dz/2, len/2, len/2+dz/2])
    true
end

@test_throws String mri_objects(:bad)

@test mri_objects_trap_test() #isa Plots.Plot
#prompt()

@test mri_objects_test2() #isa Plots.Plot
#prompt()

@test mri_objects_test_case1() #isa Plots.Plot
#prompt()

@test mri_objects_test3() #isa Plots.Plot

@test mri_objects_test_case4() #isa Plots.Plot
#prompt()
