import LinearAlgebra: lu

"""
    tripcolor(x, y, z, t; cmin=nothing, cmax=nothing)

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
    _,_,px,py = GR.inqdspsize()
    c = zeros(Int,px,py)

    xmin,xmax = extrema(x)
    ymin,ymax = extrema(y)
    isnothing(zmin) && (zmin = minimum(z))
    isnothing(zmax) && (zmax = maximum(z))
    w = xmax-xmin
    h = ymax-ymin

    nt = size(t,2)
    for it=1:nt
        xt = x[t[:,it]]
        yt = y[t[:,it]]
        zt = z[t[:,it]]

        T = [xt[1]-xt[3] xt[2]-xt[3] ; yt[1]-yt[3] yt[2]-yt[3]]
        T_lu = lu(T)

        xmint,xmaxt = extrema(xt)
        ymint,ymaxt = extrema(yt)

        x0 = floor(Int, (xmint-xmin)/w*px)
        x1 = ceil(Int, (xmaxt-xmin)/w*px)
        y0 = floor(Int, (ymint-ymin)/h*py)
        y1 = ceil(Int, (ymaxt-ymin)/h*py)

        # Loop over all pixels in triangle's bounding box
        for i=x0:x1-1, j=y0:y1-1
            # Cartesian coordinates of pixel (i,j) (zero-indexed)
            cx = xmin + i/px*w
            cy = ymin + j/py*h
            # Convert from Cartesian to barycentric coordinates
            λ1,λ2 = T_lu\[cx-xt[3], cy-yt[3]]
            λ3 = 1 - λ1 - λ2
            λ = [λ1,λ2,λ3]
            tol = 1e-12
            valid = all((λ .>= -tol) .& (λ .<= 1+tol))
            # If pixel is outside of triangle, move on to next pixel
            !valid  && continue
            zval = zt'*[λ1,λ2,λ3]
            c[i+1,py-j] = getcolorind(zval,zmin,zmax)
        end
    end
    GR.cellarray(0,1,0,1,Int(px),Int(py),c)
end
