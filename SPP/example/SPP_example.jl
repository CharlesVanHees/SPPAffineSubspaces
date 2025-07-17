include("../src/Graph.jl")
include("../src/SPP.jl")

# Create a graph with 5 vertices
G = emptyDirectedGraph(Int64)
addVertex!(G, 1)
addVertex!(G, 2)
addVertex!(G, 3)
addVertex!(G, 4)
addVertex!(G, 5)
addEdge!(G, 1, 2, 3)
addEdge!(G, 1, 4, 3)
addEdge!(G, 1, 5, 2)
addEdge!(G, 4, 5, 1)
addEdge!(G, 4, 3, 4)
addEdge!(G, 2, 3, 1)
addEdge!(G, 5, 3, 1)
addEdge!(G, 5, 1, 1)

# In this graph, the shortest path between the vertices 1 and 3 is of length 3, going through node 5.
solveSPP(G, 1, 3)