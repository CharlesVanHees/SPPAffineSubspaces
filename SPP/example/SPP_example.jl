include("../src/Graph.jl")
include("../src/SPP.jl")

# Create a graph with 5 vertices
G = emptyDirectedGraph(String)
addVertex!(G, "Apple")
addVertex!(G, "Banana")
addVertex!(G, "Peer")
addVertex!(G, "Watermelon")
addVertex!(G, "Strawberry")
addEdge!(G, "Apple", "Banana", 3)
addEdge!(G, "Apple", "Watermelon", 3)
addEdge!(G, "Apple", "Strawberry", 2)
addEdge!(G, "Watermelon", "Strawberry", 1)
addEdge!(G, "Watermelon", "Peer", 4)
addEdge!(G, "Banana", "Peer", 1)
addEdge!(G, "Strawberry", "Peer", 1)
addEdge!(G, "Strawberry", "Apple", 1)

# In this graph, the shortest path between the vertices 1 and 3 is of length 3, going through node 5.
solveSPP(G, "Apple", "Peer")