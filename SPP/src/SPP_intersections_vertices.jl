using Gurobi, JuMP

include("AffineSubspace.jl")
include("Graph.jl")
include("Problem.jl")

function problemToGraphV(P::Problem)
    I = zeros((size(P.s,1), size(P.s,1))); for i in 1:size(P.s,1) I[i,i] = 1 end # Identity matrix

    G = emptyDirectedGraph(Tuple{AffineSubspace, AffineSubspace})

    s = AffineSubspace(I, -P.s, 0.0)
    addVertex!(G, (s,s)) # The source is the first vertex
    t = AffineSubspace(I, -P.t, 0.0)
    addVertex!(G, (t,t)) # The target is the second vertex
    for i in 1:P.M, j in P.intersections[i]
        if i < j
            addVertex!(G, (P.affSubspaces[i], P.affSubspaces[j]))
            P.containSource[i] ? addEdge!(G, (s,s), (P.affSubspaces[i], P.affSubspaces[j]), P.affSubspaces[i].β) : nothing
            P.containSource[j] ? addEdge!(G, (s,s), (P.affSubspaces[i], P.affSubspaces[j]), P.affSubspaces[j].β) : nothing
            P.containTarget[i] ? addEdge!(G, (P.affSubspaces[i], P.affSubspaces[j]), (t,t), P.affSubspaces[i].β) : nothing
            P.containTarget[j] ? addEdge!(G, (P.affSubspaces[i], P.affSubspaces[j]), (t,t), P.affSubspaces[j].β) : nothing
        end
    end

    for i in 3:G.V, j in i+1:G.V
        if G.Vertices[i][1] in G.Vertices[j]
            addEdge!(G, G.Vertices[i], G.Vertices[j], G.Vertices[i][1].β)
            addEdge!(G, G.Vertices[j], G.Vertices[i], G.Vertices[i][1].β)
        end
        if G.Vertices[i][2] in G.Vertices[j]
            addEdge!(G, G.Vertices[i], G.Vertices[j], G.Vertices[i][2].β)
            addEdge!(G, G.Vertices[j], G.Vertices[i], G.Vertices[i][2].β);
        end
    end
    return G
end

function SPPGraphV(G::DirectedGraph, s::Vector{T}, t::Vector{T}; Optimizer::Module = Gurobi, R = 10000000000, verbose::Bool=true) where {T <: Real}
    n = size(s, 1)

    model = Model(Optimizer.Optimizer); set_silent(model)

    @variable(model, y_e[1:G.V, 1:G.V], Bin)
    @variables(model, begin
        x_in[ 1:G.V, 1:G.V, 1:n] # edge from i to j, from the point of view of j
        x_out[1:G.V, 1:G.V, 1:n] # edge from i to j, from the point of view of i
    end)
    for i in 1:G.V, j in 1:G.V
        if isnothing(G.Adj[i][j])      # no edge from i to j
            fix(y_e[i,j], 0)
            for k in 1:n
                fix(x_in[ i,j,k],0)
                fix(x_out[i,j,k],0)
            end
        end
    end

    @objective(model, Min, sum(G.Adj[i][j] ≠ nothing ? G.Adj[i][j] * sqrt(sum((x_out[i,j,:] .- x_in[i,j,:]).^2)) : 0 for i in 1:G.V, j in 1:G.V))

    # Flow conservation constraint 1
    @constraint(model, [i in 1:G.V], (sum(G.Adj[j][i] ≠ nothing ? y_e[j,i] : 0 for j in 1:G.V) + (i==1)) == (sum(G.Adj[i][j] ≠ nothing ? y_e[i,j] : 0 for j in 1:G.V) + (i == 2)))
    # Degree constraint
    @constraint(model, [i in 1:G.V], (sum(G.Adj[j][i] ≠ nothing ? y_e[j,i] : 0 for j in 1:G.V) + (i==1)) <= 1)

    # Flow conservation constraint 2
    @constraint(model, [i in 3:G.V], sum(G.Adj[j][i] ≠ nothing ? x_in[j,i,:] : zeros(n) for j in 1:G.V) .== sum(G.Adj[i][j] ≠ nothing ? x_out[i,j,:] : zeros(n) for j in 1:G.V))

    # Constraints of domain
    @constraints(model, begin
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[i][1].A * x_out[i,j,:] + G.Vertices[i][1].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[i][2].A * x_out[i,j,:] + G.Vertices[i][2].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[j][1].A * x_in[i,j,:] + G.Vertices[j][1].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[j][2].A * x_in[i,j,:] + G.Vertices[j][2].b * y_e[i,j]) .== 0
    end)

    # Ball
    @constraints(model, begin
        [i in 1:G.V, j in 1:G.V, k in 1:n; G.Adj[i][j] ≠ nothing], x_out[i,j,k] <=  R*y_e[i,j]
        [i in 1:G.V, j in 1:G.V, k in 1:n; G.Adj[i][j] ≠ nothing], x_out[i,j,k] >= -R*y_e[i,j]
        [i in 1:G.V, j in 1:G.V, k in 1:n; G.Adj[i][j] ≠ nothing], x_in[ i,j,k] <=  R*y_e[i,j]
        [i in 1:G.V, j in 1:G.V, k in 1:n; G.Adj[i][j] ≠ nothing], x_in[ i,j,k] >= -R*y_e[i,j]
    end)

    optimize!(model); println()
    if is_solved_and_feasible(model)
        println("A solution of optimal cost $(objective_value(model)) has been found!")
    else println("There is no path between your source and target staying on the given affine subspaces."); return
    end

    if verbose
        println("The optimal path is the following:")
        println("Source s = $(s)\n")
        i = 1
        while i != 2
            j = findfirst(Vector(value.(y_e[i,:])) .== 1)
            println("$(Vector(value.(x_out[i,j,:]))) --> $(Vector(value.(x_in[i,j,:])))")
            println("The cost is $(G.Adj[i][j] * sqrt(sum((Vector(value.(x_in[i,j,:])) - Vector(value.(x_out[i,j,:]))).^2)))\n")
            i = j
        end
        println("Target t = $(t)")
    end
end