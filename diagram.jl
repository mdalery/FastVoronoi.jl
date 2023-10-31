struct BoxedVoronoi{D, P<:Real, V<:Integer}
    points::SVector{NB_POINTS, SVector{D, P}} where NB_POINTS
    cells::Vector{Polytope{D, P, V}}
    radii::Vector{P}

    function BoxedVoronoi{V}( points::SVector{NB_POINTS, SVector{D, P}}, bounds::SVector{D, SVector{2, P}}, nb_neighbors = 10, leafsize = 10 ) where {NB_POINTS, D, P<:Real, V<:Integer}
        tree = KDTree(points, leafsize = leafsize)

        box_halfspaces = [ zeros(P, D+1) for _ in 1:2*D ]
        for k in 1:D
            lower, upper = bounds[k]

            box_halfspaces[2*k-1][k] = - one(P)
            box_halfspaces[2*k-1][end] = lower

            box_halfspaces[2*k][k] = one(P)
            box_halfspaces[2*k][end] = -upper
        end
        box_halfspaces = SVector{D+1}.(box_halfspaces)

        base_element = SVector{D}( V(2*k-1) for k in 1:D )
        box_vertices = vec([ base_element + element for element in collect.(Iterators.product([ [ zero(V), one(V) ] for _ in 1:D ]...)) ])

        cells = [ Polytope(box_halfspaces, box_vertices) for _ in 1:NB_POINTS ]
        radii = [ zero(P) for _ in 1:NB_POINTS ]
        for i in 1:NB_POINTS
            point = points[i]

            last_nearest_indices, _ = knn(tree, point, nb_neighbors, true, ind -> ind == i)
            all_nearest_indices = last_nearest_indices
            sq_radius = get_sqradius(cells[i], point)

            again = true
            while true
                for index in last_nearest_indices
                    nearest = points[index]
                    sq_distance = sum((point - nearest).^2)
                    if sq_distance > 4. * sq_radius
                        again = false
                        break
                    end

                    halfspace = vcat(2 * one(P) * (nearest - point), sum(point.^2 - nearest.^2))
                    clipped = clip!(cells[i], halfspace)
                    if clipped
                        sq_radius = get_sqradius(cells[i], point)
                    end
                end
                if !again
                    break
                end
                last_nearest_indices, _ = knn(tree, point, nb_neighbors, true, ind -> ind in [ i; all_nearest_indices ])
                all_nearest_indices = [ last_nearest_indices; all_nearest_indices ]
                if length(all_nearest_indices) == NB_POINTS - 1
                    break
                end
            end

            radii[i] = sqrt(sq_radius)
        end

        return new{D, P, V}(points, cells, radii)
    end
end
