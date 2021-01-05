import GR
using Images, Plots, Triplot

include("common.jl")

function to_rgba(x::UInt32)
    a = ((x & 0xff000000)>>24)/255
    b = ((x & 0x00ff0000)>>16)/255
    g = ((x & 0x0000ff00)>>8)/255
    r = (x & 0x000000ff)/255
    RGBA(r, g, b, a)
end

x,y,t = cartmesh(100)
z = peaks_unit.(x,y)

GR.setcolormap(1)
img = to_rgba.(Triplot.rasterize(x,y,z,t)')
plot([extrema(x)...],[extrema(y)...],img,yflip=false)
