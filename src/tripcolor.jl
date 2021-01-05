"""
    tripcolor(x, y, z, t; zmin=nothing, zmax=nothing)

Draw a pseudocolor plot of unstructured triangular data. `x`, `y`, and `z` are
one-dimensional input arrays of the same length (the number of points). `t` is an integer
array of size `(3,nt)`, where `nt` is the number of triangles. The coordinates of the `i`th
vertex of the `j`th triangle are given by `(x[t[i,j]], y[t[i,j]])`. Function values at
vertices are specified by the `z` array

GR's current colormap (set using `GR.setcolormap`) is used. The colormap bounds can be set
by specifying values for `zmin` and `zmax`. If these parameters are set to `nothing`, then
the minimum and maximum values of the `z` array will be used.
"""
function tripcolor(x, y, z, t; zmin=nothing, zmax=nothing)
    c = rasterize(x, y, z, t; zmin, zmax)
    GR.drawimage(0,1,0,1,size(c,1),size(c,2),c)
end

function rasterize(x, y, z, t; px=nothing, py=nothing, zmin=nothing, zmax=nothing)
    _,_,px_dsp,py_dsp = GR.inqdspsize()
    isnothing(px) && (px = px_dsp)
    isnothing(py) && (py = py_dsp)

    c = zeros(UInt32,px,py)

    # Precompute encoded RGBA values of colormap values
    cmap = zeros(UInt32,256)
    for i=0:255
        # Pixel is encoded as 0xAABBGGRR, GR.inqcolor returns 0x00BBGGRR, so set opaque
        cmap[begin+i] = 0xff000000 | GR.inqcolor(1000+i)
    end

    xmin,xmax = extrema(x)
    ymin,ymax = extrema(y)
    w = xmax-xmin
    h = ymax-ymin
    isnothing(zmin) && (zmin = minimum(z))
    isnothing(zmax) && (zmax = maximum(z))

    nt = size(t,2)
    for it=1:nt
        xt = @views x[t[:,it]]
        yt = @views y[t[:,it]]
        zt = @views z[t[:,it]]

        detT = (yt[2]-yt[3])*(xt[1]-xt[3]) + (xt[3]-xt[2])*(yt[1]-yt[3])

        xmint,xmaxt = extrema(xt)
        ymint,ymaxt = extrema(yt)

        x0 = floor(Int,(xmint-xmin)/w*px)
        x1 = ceil(Int,(xmaxt-xmin)/w*px)
        y0 = floor(Int,(ymint-ymin)/h*py)
        y1 = ceil(Int,(ymaxt-ymin)/h*py)

        # Loop over all pixels in triangle's bounding box
        for i=x0:x1-1, j=y0:y1-1
            # Cartesian coordinates of pixel (i,j) (zero-indexed)
            cx = xmin + i/px*w
            cy = ymin + j/py*h
            # Convert from Cartesian to barycentric coordinates
            λ1 = ((yt[2] -yt[3])*(cx-xt[3]) + (xt[3]-xt[2])*(cy-yt[3]))/detT
            λ2 = ((yt[3] -yt[1])*(cx-xt[3]) + (xt[1]-xt[3])*(cy-yt[3]))/detT
            λ3 = 1 - λ1 - λ2
            tol = 1e-12
            # If pixel is outside of triangle, move on to next pixel
            (λ1>-tol && λ2<1+tol && λ2>-tol && λ2<1+tol && λ3>-tol && λ3<1+tol) || continue
            zval = zt[1]*λ1 + zt[2]*λ2 + zt[3]*λ3
            zval = clamp(zval,zmin,zmax)
            c[i+1,py-j] = cmap[begin+getcolorind(zval,zmin,zmax)-1000]
        end
    end
    c
end
