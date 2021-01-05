import GR
using Images, Plots, Triplot

include("common.jl")

x,y,t = cartmesh(100)
GR.setcolormap(1)
img = reinterpret(RGB24,Triplot.rasterize(x,y,peaks_unit.(x,y),t))
plot([0,1],[0,1],img,yflip=false)
