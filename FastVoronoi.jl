module FastVoronoi

using StaticArrays,
      NearestNeighbors,
      LinearAlgebra,
      Plots

include("polytope.jl")

include("diagram.jl")
export BoxedVoronoi

include("2dplots.jl")
export plot

end
