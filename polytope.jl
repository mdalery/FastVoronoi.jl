mutable struct Polytope{D, H<:Real, V<:Integer}
    halfspaces::Vector{SVector{N, H}} where N
    vertices::Vector{SVector{D, V}}

    function Polytope( halfspaces::Vector{SVector{N, H}}, vertices::Vector{SVector{D, V}} ) where {N, D, H<:Real, V<:Integer}
        @assert N == D+1
        return new{D, H, V}(copy(halfspaces), copy(vertices))
    end
end

function dimension( p::Polytope{D} ) where D
    return D
end

function intersection( p::Polytope{D, <:Real, V}, vertex::SVector{D, V} ) where {D, V<:Integer}
    mat = hcat(p.halfspaces[vertex]...)'

    return mat[:, 1:end-1] \ (-mat[:, end])
end

function clip!( p::Polytope{D, H, V}, halfspace::SVector{N, H} ) where {N, D, H<:Real, V<:Integer}
    @assert N == D+1
    removed_vertices = SVector{D, V}[]
    for vertex in p.vertices
        if halfspace ⋅ vcat(intersection(p, vertex), one(H)) > zero(H)
            push!(removed_vertices, vertex)
        end
    end
    if length(removed_vertices) == 0
        return false
    end

    p.vertices = [ vertex for vertex in p.vertices if vertex ∉ removed_vertices ]
    push!(p.halfspaces, halfspace)
    index_halfspace = length(p.halfspaces)

    add_vertices = SVector{D, V}[]
    for r in removed_vertices
        for k in p.vertices
            edge = intersect(r, k)
            if length(edge) == D-1
                push!(add_vertices, [ edge; index_halfspace ])
            end
        end
    end
    append!(p.vertices, add_vertices)

    return true
end

function get_sqradius( p::Polytope{D, H}, center::SVector{D, H} ) where {D, H<:Real}
    sq_radius = zero(H)

    for vertex in p.vertices
        vertex_coordinates = intersection(p, vertex)
        sq_distance = sum((vertex_coordinates - center).^2)
        if sq_distance > sq_radius
            sq_radius = sq_distance
            furthest_vertex = vertex
        end
    end

    return sq_radius
end
