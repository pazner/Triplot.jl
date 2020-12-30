const XY = Tuple{Float64, Float64}
const Contour = Vector{XY}

"""
    tricontour(x, y, z, t, levels)

Draw contour lines of unstructured triangular data. `x`, `y`, and `z` are one-dimensional
input arrays of the same length (the number of points). `t` is an integer array of size
`(3,nt)`, where `nt` is the number of triangles. The coordinates of the `i`th vertex of the
`j`th triangle are given by `(x[t[i,j]], y[t[i,j]])`. Function values at vertices are
specified by the `z` array, and contours levels are defined by `levels.`

If `levels` is an integer, then equally spaced levels between the extrema of `z` will be
generated. Otherwise, `levels` is a vector of increasing values, each value corresponding
to one contour line.

The contour lines will be colored using the current colormap. Contour lines of a single
color lines can be obtained by first setting a completely uniform colormap, for example

    GR.setcolormapfromrgb([0.0,0.0], [0.0,0.0], [0.0,0.0])
"""
function tricontour(x, y, z, t, levels)
    m = TriMesh(x, y, t)
    if levels isa Integer
        levels = range(minimum(z), maximum(z), length=levels)
    end
    tricontour(m, z, levels)
end

function tricontour(m::TriMesh, z, levels)
    @assert issorted(levels)
    lmin,lmax = extrema(levels)
    for level=levels
        contours = generate_unfilled_contour(m, z, level)
        GR.setlinecolorind(getcolorind(level,lmin,lmax))
        for c=contours
            x = first.(c)
            y = last.(c)
            GR.polyline(x,y)
        end
    end
end

function getcolorind(level,lmin,lmax)
    cfrac = (lmin == lmax) ? 0 : (level-lmin)/(lmax-lmin)
    round(Int,1000 + 255*cfrac)
end

function generate_unfilled_contour(m::TriMesh, z, level)
    reset_visited(m)
    contours = Vector{Contour}()
    add_bdr_lines!(contours, m, z, level)
    add_interior_lines!(contours, m, z, level; on_upper=false, filled=false)
    contours
end

function interp(m::TriMesh, p1, p2, z, level)
    t = (z[p2] - level)/(z[p2] - z[p1])
    x = t*m.x[p1] + (1.0-t)*m.x[p2]
    y = t*m.y[p1] + (1.0-t)*m.y[p2]
    x,y
end

function interp(m::TriMesh, edge::TriEdge, z, level)
    p1 = m.t[edge.ie,edge.it]
    p2 = m.t[mod1(edge.ie+1,3),edge.it]
    interp(m, p1, p2, z, level)
end

function get_exit_edge(m::TriMesh, it, z, level; on_upper)
    config = (z[m.t[1,it]] >= level) |
             (z[m.t[2,it]] >= level) << 1 |
             (z[m.t[3,it]] >= level) << 2
    on_upper && (config = 7-config)
    if config == 0 return 0
    elseif config == 1 return 3
    elseif config == 2 return 1
    elseif config == 3 return 3
    elseif config == 4 return 2
    elseif config == 5 return 2
    elseif config == 6 return 1
    elseif config == 7 return 0
    else error("Invalid config value")
    end
end

function follow_interior!(contour::Contour, m::TriMesh, edge::TriEdge, z, level;
                          end_on_bdr, on_upper=false)
    push!(contour, interp(m, edge, z, level))
    it = edge.it
    while true
        vi = on_upper ? it + size(m.t,2) : it
        !end_on_bdr && m.visited[vi] && return edge
        m.visited[vi] = true
        ie = get_exit_edge(m, it, z, level; on_upper)
        @assert ie > 0
        edge = TriEdge(ie,it)
        push!(contour, interp(m, edge, z, level))
        # Move to next triangle
        next_edge = get_neighbor_edge(m,edge)
        it = next_edge.it
        end_on_bdr && it == 0 && return edge
        edge = next_edge
    end
end

function follow_interior(m::TriMesh, edge::TriEdge, z, level; end_on_bdr, on_upper=false)
    contour = Contour()
    edge = follow_interior!(contour, m, edge, z, level; end_on_bdr, on_upper)
    contour, edge
end

function add_bdr_lines!(contours, m::TriMesh, z, level)
    for bdr=m.bdr_loops
        nedges = length(bdr)
        end_above = false
        for i=1:nedges
            edge = bdr[i]
            if i == 1
                start_above = z[m.t[edge.ie,edge.it]] >= level
            else
                start_above = end_above
            end
            end_above = z[m.t[mod1(edge.ie+1,3),edge.it]] >= level
            if start_above && !end_above
                contour,_ = follow_interior(m, edge, z, level; end_on_bdr=true)
                push!(contours, contour)
            end
        end
    end
end

function add_interior_lines!(contours, m::TriMesh, z, level; on_upper, filled)
    nt = size(m.t,2)
    for it=1:nt
        vi = on_upper ? it + size(m.t,2) : it
        m.visited[vi] && continue
        m.visited[vi] = true
        ie = get_exit_edge(m, it, z, level; on_upper)
        ie == 0 && continue
        # Found a new outgoing contour line, start from neighbor
        edge = get_neighbor_edge(m, TriEdge(ie,it))
        contour,_ = follow_interior(m, edge, z, level; on_upper, end_on_bdr=false)
        push!(contour, first(contour)) # Make closed loop
        push!(contours, contour)
    end
end
