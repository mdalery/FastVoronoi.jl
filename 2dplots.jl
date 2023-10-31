function plot( vor::BoxedVoronoi )
    scatter(first.(vor.points), last.(vor.points), markersize = 1, color = "blue", label = "")

    for i in 1:length(vor.cells)
        vertex = vor.cells[i].vertices[1]
        vertices = [ vertex ]
        coords_vertices = [ FastVoronoi.intersection(vor.cells[i], vertex) ]

        for k in 1:length(vor.cells[i].vertices)-1
            vertex = vor.cells[i].vertices[findfirst(x -> x ∉ vertices && length(intersect(x, vertex)) == 1, vor.cells[i].vertices)]
            push!(coords_vertices, FastVoronoi.intersection(vor.cells[i], vertex))
            push!(vertices, vertex)
        end
        push!(coords_vertices, coords_vertices[1])
        display(plot!(first.(coords_vertices), last.(coords_vertices), color = "red", label = ""))
    end
end
