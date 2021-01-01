import GR
using Tricontour

include("common.jl")

function plot_pcolor(fname, n, f)
    x,y,t = cartmesh(n)
    z = f.(x,y)

    GR.beginprint(fname)
    GR.setcolormap(1)
    tripcolor(x,y,z,t)
    GR.endprint()
end

plot_pcolor(joinpath(@__DIR__, "diag_pcolor.png"), 5, (x,y) -> x + y)
plot_pcolor(joinpath(@__DIR__, "peaks_pcolor.png"), 100, peaks_unit)
plot_pcolor(joinpath(@__DIR__, "sine_pcolor.png"), 100, (x,y) -> sin(2*pi*x)*cos(2*pi*y))
