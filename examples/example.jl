import GR
using Tricontour

meshgrid(x,y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

function cartmesh(n)
    nv = n^2
    x1d = range(0,1,length=n)
    x,y = meshgrid(x1d,x1d)
    nt = 2*(n-1)^2
    t = zeros(Int, 3, nt)
    it = 1
    iv = 1
    for iy=1:n-1
        for ix=1:n-1
            t[:,it] = [iv,iv+1,iv+n]
            t[:,it+1] = [iv+1,iv+n+1,iv+n]
            it += 2
            iv += 1
        end
        iv += 1
    end
    x,y,t
end

peaks(x,y) = (3*(1-x)^2*exp(-(x^2) - (y+1)^2)
    - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2)
    - 1/3*exp(-(x+1)^2 - y^2))
peaks_unit(x,y) = peaks(2*pi*(x-0.5), 2*pi*(y-0.5))

function plot_contours(fname, n, f; levels=32)
    x,y,t = cartmesh(n)
    z = f.(x,y)

    GR.beginprint(fname)
    GR.setcolormap(1)
    GR.setlinewidth(1.0)
    tricontourf(x,y,z,t,levels)
    GR.setcolormapfromrgb([0.0,0.0],[0.0,0.0],[0.0,0.0])
    tricontour(x,y,z,t,levels)
    GR.endprint()
end

plot_contours(joinpath(@__DIR__, "diag.png"), 5, (x,y) -> x + y)
plot_contours(joinpath(@__DIR__, "peaks.png"), 100, peaks_unit)
plot_contours(joinpath(@__DIR__, "sine.png"), 100, (x,y) -> sin(2*pi*x)*cos(2*pi*y))
