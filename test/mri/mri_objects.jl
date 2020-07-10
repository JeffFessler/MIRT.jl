# mri_objects.jl

using MIRT: mri_objects

using MIRT: jim, image_geom, prompt, max_percent_diff
import MIRT: mri_objects_trap
using FFTW: fft, fftshift, ifftshift
using Plots
using Test: @test


# 2D test
#NOTE: Dirac will not plot very well.
#Just make sure it prints a gray image, and that's fine. FFT will also be screwed up.
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
        isinteractive() && (@show max_percent_diff(f2,s2))
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
    xt = mri_objects(:case1,unit = :cm).image(ig.xg, ig.yg)
    jim(xt, title = "case1cm")
end


function mri_objects_test_case4()
    ig = image_geom(nx = 2^5, ny = 2^5 + 2,nz = 2^5, dx = 1,dz = .5, offsets = :dsp)
    xt = mri_objects(:case4,unit = :cm).image(ig.xg, ig.yg, ig.zg)
    jim(xt, title = "case4cm")
end


# 3D test
# NOTE: Dirac will not plot very well.
# Just make sure it prints a gray image, and that's fine. FFT will also be screwed up.
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
        isinteractive() && (@show max_percent_diff(f3,s3))
        pl[i,1] = jim(i3, title = "$(test[1])")
        pl[i,2] = jim(abs.(s3), title = "fft")
        pl[i,3] = jim(abs.(f3), title = "kspace")
        plot(pl[i,:]...)
        prompt()
    end
    return plot(pl...)
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


@test mri_objects_trap_test() isa Plots.Plot
prompt()

@test mri_objects_test2() isa Plots.Plot
prompt()

@test mri_objects_test_case1() isa Plots.Plot
prompt()

@test mri_objects_test3() isa Plots.Plot

@test mri_objects_test_case4() isa Plots.Plot
prompt()
