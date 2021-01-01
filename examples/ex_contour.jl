import GR
using Tricontour

include("common.jl")

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
