using Gurobi, JuMP

include("AffineSubspace.jl")
include("Graph.jl")
include("Problem.jl")

function problemToGraphE(P::Problem)
    I = zeros((size(P.s,1), size(P.s,1))); for i in 1:size(P.s,1) I[i,i] = 1 end # Identity matrix

    G = emptyDirectedGraph(AffineSubspace)

    s = AffineSubspace(I, -P.s, 0.0) # The source can be seen as an affine subspace of equation x - s = 0.
    addVertex!(G, s) # The first vertex is the source
    for AS in P.affSubspaces addVertex!(G, AS) end
    t = AffineSubspace(I, -P.t, 0.0)
    addVertex!(G, t) # The last vertex is the target

    for i in findall(P.containSource) addEdge!(G, s, P.affSubspaces[i], 0.0) end
    for i in 1:P.M
        for e in P.intersections[i] addEdge!(G, P.affSubspaces[i], P.affSubspaces[e], P.affSubspaces[i].β) end
    end
    for i in findall(P.containTarget) addEdge!(G, P.affSubspaces[i], t, P.affSubspaces[i].β) end
    return G
end

function SPPGraphE(G::DirectedGraph, s::Vector{T}, t::Vector{T}, Optimizer = Gurobi, verbose::Bool=true) where {T <: Real}
    n = size(G.Vertices[1].A, 2)

    model = Model(Optimizer.Optimizer); set_silent(model)

    @variable(model, y_e[1:G.V, 1:G.V], Bin)
    @variables(model, begin
        q_0_in[ 1:G.V, 1:G.V, 1:n] # edge from i to j, from the point of view of j
        q_0_out[1:G.V, 1:G.V, 1:n] # edge from i to j, from the point of view of i
        q_1_in[ 1:G.V, 1:G.V, 1:n]
        q_1_out[1:G.V, 1:G.V, 1:n]
    end)
    for i in 1:G.V, j in 1:G.V
        if isnothing(G.Adj[i][j])       # no edge between i and j
            fix(y_e[i,j], 0)
            for k in 1:n
                fix(q_0_in[ i,j,k], 0)
                fix(q_0_out[i,j,k], 0)
                fix(q_1_in[ i,j,k], 0)
                fix(q_1_out[i,j,k], 0)
            end
        end
    end

    @objective(model, Min, sum(G.Adj[i][j] ≠ nothing ? G.Adj[i][j] * sqrt(sum((q_1_out[i,j,:] .- q_0_out[i,j,:]).^2)) : 0 for i in 2:G.V, j in 2:G.V))

    # Flow conservation constraint 1
    @constraint(model, [i in 1:G.V], (sum(G.Adj[j][i] ≠ nothing ? y_e[j,i] : 0 for j in 1:G.V) + (i==1)) == (sum(G.Adj[i][j] ≠ nothing ? y_e[i,j] : 0 for j in 1:G.V) + (i==G.V)))
    # Degree constraint
    @constraint(model, [i in 1:G.V], (sum(G.Adj[j][i] ≠ nothing ? y_e[j,i] : 0 for j in 1:G.V) + (i==1)) <= 1)

    # Flow conservation constraint 2
    @constraint(model, [i in 2:G.V-1], sum(G.Adj[j][i] ≠ nothing ? q_0_in[j,i,:] : zeros(n) for j in 1:G.V) .== sum(G.Adj[i][j] ≠ nothing ? q_0_out[i,j,:] : zeros(n) for j in 1:G.V))
    @constraint(model, [i in 2:G.V-1], sum(G.Adj[j][i] ≠ nothing ? q_1_in[j,i,:] : zeros(n) for j in 1:G.V) .== sum(G.Adj[i][j] ≠ nothing ? q_1_out[i,j,:] : zeros(n) for j in 1:G.V))

    # Constraints of domain
    @constraints(model, begin
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[i].A * q_0_out[i,j,:] + G.Vertices[i].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[i].A * q_1_out[i,j,:] + G.Vertices[i].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[j].A * q_0_in[ i,j,:] + G.Vertices[j].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[j].A * q_1_in[ i,j,:] + G.Vertices[j].b * y_e[i,j]) .== 0
    end)

    # Continuity constraint
    @constraint(model, [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_1_out[i,j,:] .== q_0_in[i,j,:])

    # Ball
    @constraints(model, begin
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_0_out[i,j,:] .<= 10*ones(n)*y_e[i,j]
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_0_out[i,j,:] .>= -10*ones(n)*y_e[i,j]
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_1_out[i,j,:] .<= 10*ones(n)*y_e[i,j]
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_1_out[i,j,:] .>= -10*ones(n)*y_e[i,j]
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_0_in[i,j,:] .<= 10*ones(n)*y_e[i,j]
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_0_in[i,j,:] .>= -10*ones(n)*y_e[i,j]
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_1_in[i,j,:] .<= 10*ones(n)*y_e[i,j]
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], q_1_in[i,j,:] .>= -10*ones(n)*y_e[i,j]
    end)

    optimize!(model); println()
    if is_solved_and_feasible(model)
        println("A solution of optimal cost $(objective_value(model)) has been found!")
    else println("There is no path between your source and target staying on the given affine subspaces."); return
    end

    if verbose
        println("The optimal path is the following:")
        i = findfirst(Vector(value.(y_e[1,:])) .== 1)
        println("Source s = $(s) on AS $(i - 1)")
        while i != G.V
            j = findfirst(Vector(value.(y_e[i,:])) .== 1)
            print("$(Vector(value.(q_0_out[i,j,:]))) --> $(Vector(value.(q_1_out[i,j,:]))) to reach ")
            j == G.V ? println("target t") : println("AS $(j-1)")
            println("The cost is $(G.Adj[i][j] * sqrt(sum((Vector(value.(q_0_out[i,j,:])) - Vector(value.(q_1_out[i,j,:]))).^2)))")
            i = j
        end
        println("Target t = $(t)")
    end
end