using Printf, Random

mutable struct DirectedGraph{T}
    V::UInt64 # Number of Vertices
    E::UInt64 # Number of Edges

    Vertices::Vector{T}     # Vertices
    Adj::Vector{Vector{Union{Nothing, Number}}}  # Adjency Matrix of the graph - this structure assume a dense graph
end

emptyDirectedGraph(T::Union{DataType, UnionAll}) = DirectedGraph(UInt64(0), UInt64(0), Vector{T}(undef, 0), Vector{Vector{Union{Nothing, Number}}}(undef, 0))

function addVertex!(G::DirectedGraph, value::T) where {T}
    @assert T <: eltype(G.Vertices)
    if (value in G.Vertices) return end
    for i in eachindex(G.Adj) push!(G.Adj[i], nothing) end
    G.V += 1
    push!(G.Vertices, value)
    push!(G.Adj, Vector{eltype(eltype(G.Adj))}(nothing, G.V))
end

function addEdge!(G::DirectedGraph, i::T, j::T, w::Number) where {T}
    @assert T <: eltype(G.Vertices)
    i, j = findfirst(==(i), G.Vertices), findfirst(==(j), G.Vertices)
    @assert i ≠ nothing && j ≠ nothing
    if (G.Adj[i][j] ≠ nothing)
        G.Adj[i][j] = min(G.Adj[i][j], w) # We take the minimum between the two weights
    else
        G.E += 1; G.Adj[i][j] = w
    end
end

function printGraph(G::DirectedGraph)
    println("Number of vertices: $(G.V)")
    println("Number of edges: $(G.E)\n")

    println("Vertices list:")
    print("[ "); for i in eachindex(G.Vertices) print("$(G.Vertices[i]) ") end; println("]")

    println("Adjacency matrix:")
    for i in eachindex(G.Adj)
        print("Vertex $(i): [\t")
        for j in G.Adj[i] print("$(j)\t") end
        println("]")
    end
end

function createRandomDirectedGraph()
    G = emptyDirectedGraph(Int64)
    for i in 1:rand(1:10) addVertex!(G, i) end
    for i in 1:G.V
        for _ in 1:rand(1:5) addEdge!(G, Int64(i), Int64(rand(1:G.V)), rand(Float64)) end
    end
    return G
end

# G = createRandomDirectedGraph()
# printGraph(G)