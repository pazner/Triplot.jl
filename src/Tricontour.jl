module Tricontour

using GR

export tricontour, tricontourf

include("trimesh.jl")
include("contour.jl")
include("contourf.jl")

end
