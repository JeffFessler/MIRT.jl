using Plots

"""
`(sino, pos, ang) = rect_sino(sg, rects;
            oversample=1, xscale=1, yscale=1)`

Create sinogram projection of one or more rectangles.
Works for both parallel-beam geometry and for fan-beam geometry.

in
* `sg`                  sinogram geometry object from sino_geom()
* `rects`               [ne 6]
                        [centx centy widthx width y angle_degrees value]

options
* `oversample::Integer` oversampling factor for emulating "strips"
                        default 1: just 1 ray per detector element
* `xscale::Integer`     use -1 to flip in x (not recommended)
* `yscale::Integer`     use -1 to flip in y (not recommended)
* `type::Symbol`        `:fan` `:par` `:moj`

out
* `sino`                [nb na] sinogram
* `pos`                 [nb 1]  radial sample positions (or [nb na] is mojette)
* `ang`                 [na]    angular sample locations (in radians)
"""
function rect_sino(sg::MIRT_sino_geom,
    ells::AbstractArray{<:Real,1};
    oversample::Integer=1,
    xscale::Integer=1,
    yscale::Integer=1,
    type::Symbol=:fan)

    type = sg.type
    if type == :fan
        (sino, pos) = rect_sino_go(rects, sg.nb, sg.ds, sg.offset_s, sg.na,
                        sg.ar, sg.dso, sg.dod, sg.dfs, sg.source_offset,
                        xscale, yscale, oversample, 0)
    elseif type == :par
        (sino, pos) = rect_sino_go(rects, sg.nb, sg.dr, sg.offset_r, sg.na,
                        sg.ar, Inf, 1, 0, sg.source_offset, xscale, yscale, oversample, 0)
    elseif type == :moj
        (sino, pos) = rect_sino_go(rects, sg.nb, [], sg.offset_r, sg.na,
                        sg.ar, Inf, 1, 0, sg.source_offset, xscale, yscale, oversample, 0)
    else
        throw("sino geom $type not done")
    end

    return (sino, pos, ang)
end

"""
`rect_sino_go()`
"""
function rect_sino_go(rects, nb, ds, offset_s, na, ang, dso, dod, dfs,
    source_offset, xscale, yscale, oversample, mojette)

    (pos, pos2) = rect_sino_pos([], nb, ds, offset_s, oversample, mojette, ang)

    sino = rect_sino_do(rects, pos2, ang[:]', xscale, yscale, dso, dod, dfs, source_offset)

    if oversample != 0
        sino = downsample2(sino, [oversample 1])
    end
    return (sino, pos)
end

"""
`rect_sino_pos()`

determine usual and fine "radial" sample positions
"""
function rect_sino_pos(pos, nb, ds, offset_s, nover, mojette, ang)

    if isempty(pos)
        if isempty(nb)
            throw("nb required when no pos provided")
        end
        wb = (nb - 1)/2 + offset_s
    else
        if !isempty(nb)
            throw("nb ignored when pos provided")
        end
    end

    if mojette != 0 # tricky mojette radial sampling
        # trick: ray_spacing aka ds comes from dx which is in mojette
        if !isempty(ds)
            throw("ds must be empty for mojette case")
        end
        dt = abs(mojette) * max(abs(cos(ang)), abs(sin(ang)))' # [1, na]

        if !isempty(pos)
            throw("mojette requires empty 'pos'")
        end

        na = max(size(ang))
        pos_coarse = ([0:nb - 1] - wb)' * dt[:]'
        if nover > 1
            pos_fine = zeros(nover*nb, na)
            for ia in 1:na
                tmp = [-(over - 1):2:(nover - 1)]' / (2*nover) * dt[ia]
                #pos_fine[:, ia] = col(outer_sum(tmp, pos_coarse[:, is]'))
            end
        else
            pos_fine = pos_coarse
        end
    else # ordinary mojette sampling
        if isempty(pos)
            pos_coarse = ([0:nb - 1]' - wb) * ds
        else
            pos_coarse = pos[:]
            ds = pos[2] - pos[1]
            #if any(abs(diff(pos) / ds - 1) > 1e-10)
        #        throw("uniform spacing required")
        #    end
        end
        if nover > 1
            # determine fine sampling positions
            pos_fine = [-(nover - 1):2:(nover-1)]' / (2*nover) * ds
            pos_fine = outer_sum(pos_fine, pos_coarse[:]')
            pos_fine = pos_fine[:]
        else
            pos_fine = pos_coarse
        end
    end
    return (pos_coarse, pos_fine)
end

"""
`rect_sino_do()`

analytical line-integral projections rectangle
"""
function rect_sino_do(rects, pos, ang, xscale, yscale, dso, dod, dfs, source_offset)


end
