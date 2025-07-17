using Gurobi, JuMP

include("AffineSubspace.jl")
include("Graph.jl")
include("Problem.jl")

myOptimizer = Gurobi

function problemToGraphE(P::Problem)
    G = emptyDirectedGraph(AffineSubspace)
    for AS in P.affSubspaces addVertex!(G, AS) end
    for i in 1:P.M
        for e in P.intersections[i] addEdge!(G, P.affSubspaces[i], P.affSubspaces[e], P.affSubspaces[i].β) end
    end
    return G
end

function SPPGraphE(G::DirectedGraph, s::Vector{T}, t::Vector{T}, verbose::Bool=true) where {T <: Real}
    n = size(G.Vertices[1].A, 2)

    model = Model(Gurobi.Optimizer)
    # set_silent(model)

    @variable(model, y_e[1:G.V, 1:G.V], Bin)
    @variables(model, begin
        q_0_in[ 1:G.V, 1:G.V, 1:n] # edge from i to j, from the point of view of j
        q_0_out[1:G.V, 1:G.V, 1:n] # edge from i to j, from the point of view of i
        q_1_in[ 1:G.V, 1:G.V, 1:n]
        q_1_out[1:G.V, 1:G.V, 1:n]
    end)
    for i in 1:G.V, j in 1:G.V              # no edge between i and j
        if isnothing(G.Adj[i][j])
            fix(y_e[i,j], 0)
            for k in 1:n
                fix(q_0_in[ i,j,k], 0)
                fix(q_0_out[i,j,k], 0)
                fix(q_1_in[ i,j,k], 0)
                fix(q_1_out[i,j,k], 0)
            end
        end
    end

    @objective(model, Min, sum(G.Adj[i][j] ≠ nothing ? G.Adj[i][j] * sum((q_1_out[i,j,:] .- q_0_out[i,j,:]).^2) : 0 for j in 1:G.V, i in 1:G.V))

    # Flow conservation constraint
    @constraint(model, [i in 1:G.V], (sum(y_e[:,i]) + (contains(G.Vertices[i], s))) == (sum(y_e[i,:]) + (contains(G.Vertices[i], t))))
    # Degree constraint
    @constraint(model, [i in 1:G.V], (sum(y_e[:,i])) + (contains(G.Vertices[i], s)) <= 1)

    # Constraint ...
    @constraint(model, [i in 1:G.V], sum(G.Adj[i][j] ≠ nothing ? q_0_in[j,i,:] : zeros(n) for j in 1:G.V) .== sum(G.Adj[i][j] ≠ nothing ? q_0_out[i,j,:] : zeros(n) for j in 1:G.V))

    # Constraint ...
    @constraints(model, begin
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], G.Vertices[i].A * q_0_out[i,j,:] + G.Vertices[i].b * y_e[i,j] == 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], G.Vertices[i].A * q_1_out[i,j,:] + G.Vertices[i].b * y_e[i,j] == 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], G.Vertices[j].A * q_0_in[ i,j,:] + G.Vertices[j].b * y_e[i,j] == 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], G.Vertices[j].A * q_1_in[ i,j,:] + G.Vertices[j].b * y_e[i,j] == 0
    end)

    # Continuity constraint
    @constraint(model, [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_1_out[i,j,:] .== q_0_in[i,j,:])

    relax_integrality(model)
    optimize!(model)
end

function example()
    P = exampleProblem()
    G = problemToGraphE(P)
    printGraph(G)
    SPPGraphE(G, P.s, P.t)
end

example()