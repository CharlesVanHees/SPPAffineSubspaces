using Gurobi, JuMP

include("AffineSubspace.jl")
include("Graph.jl")
include("Problem.jl")

function problemToGraphV(P::Problem)
    I = zeros(Int64, (size(P.s,1), size(P.s,1))); for i in 1:size(P.s,1) I[i,i] = 1 end # Identity matrix

    G = emptyDirectedGraph(Tuple{AffineSubspace, AffineSubspace, Int64, Int64})

    s = AffineSubspace(I, -P.s, 0.)
    addVertex!(G, (s,s,0,0)) # The source is the first vertex
    t = AffineSubspace(I, -P.t, 0.)
    addVertex!(G, (t,t,0,0)) # The target is the second vertex
    for i in 1:P.M, j in P.intersections[i]
        if i < j
            addVertex!(G, (P.affSubspaces[i], P.affSubspaces[j], i, Int64(j)))
            P.containSource[i] ? addEdge!(G, (s,s,0,0), (P.affSubspaces[i], P.affSubspaces[j],i,Int64(j)), P.affSubspaces[i].β) : nothing
            P.containSource[j] ? addEdge!(G, (s,s,0,0), (P.affSubspaces[i], P.affSubspaces[j],i,Int64(j)), P.affSubspaces[j].β) : nothing
            P.containTarget[i] ? addEdge!(G, (P.affSubspaces[i], P.affSubspaces[j],i,Int64(j)), (t,t,0,0), P.affSubspaces[i].β) : nothing
            P.containTarget[j] ? addEdge!(G, (P.affSubspaces[i], P.affSubspaces[j],i,Int64(j)), (t,t,0,0), P.affSubspaces[j].β) : nothing
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
    for i in 1:P.M
        P.containSource[i] && P.containTarget[i] ? addEdge!(G, (s,s,0,0), (t,t,0,0), P.affSubspaces[i].β) : nothing
    end

    return G
end

function SPPGraphV(G::DirectedGraph, s::Vector, t::Vector; Optimizer::Module = Gurobi, R = 1000, verbose::Bool=true)
    n = size(s, 1)

    model = Model(Optimizer.Optimizer)# ; set_silent(model)
    if !verbose set_silent(model) end
    # set_optimizer_attribute(model, "NumericFocus", 3)

    @variable(model, y_e[1:G.V, 1:G.V], Bin)
    @variables(model, begin
        x_in[ 1:G.V, 1:G.V, 1:n] # edge from i to j, from the point of view of j
        x_out[1:G.V, 1:G.V, 1:n] # edge from i to j, from the point of view of i
    end)
    @variable(model, w[1:G.V, 1:G.V]) # Auxiliary variables for the second order cone
    for i in 1:G.V, j in 1:G.V
        if isnothing(G.Adj[i][j])     # no edge from i to j
            fix(y_e[i,j], 0)
            for k in 1:n
                fix(x_in[ i,j,k],0)
                fix(x_out[i,j,k],0)
            end
            fix(w[i,j], 0)
        end
    end

    @objective(model, Min, sum(G.Adj[i][j] ≠ nothing ? G.Adj[i][j] * w[i,j] : 0 for i in 1:G.V, j in 1:G.V))

    # Second order cone constraint
    @constraint(model, [i in 1:G.V, j in 1:G.V], [w[i,j]; (x_out[i,j,:] .- x_in[i,j,:])] in SecondOrderCone())

    # Flow conservation constraint 1
    @constraint(model, [i in 1:G.V], (sum(G.Adj[j][i] ≠ nothing ? y_e[j,i] : 0 for j in 1:G.V) + (i==1)) == (sum(G.Adj[i][j] ≠ nothing ? y_e[i,j] : 0 for j in 1:G.V) + (i==2)))
    # Degree constraint
    @constraint(model, [i in 1:G.V], (sum(G.Adj[j][i] ≠ nothing ? y_e[j,i] : 0 for j in 1:G.V) + (i==1)) <= 1)

    # Flow conservation constraint 2
    @constraint(model, [i in 3:G.V], sum(G.Adj[j][i] ≠ nothing ? x_in[j,i,:] : zeros(n) for j in 1:G.V) .== sum(G.Adj[i][j] ≠ nothing ? x_out[i,j,:] : zeros(n) for j in 1:G.V))

    # Constraints of domain
    @constraints(model, begin
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[i][1].A * x_out[i,j,:] + G.Vertices[i][1].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[i][2].A * x_out[i,j,:] + G.Vertices[i][2].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[j][1].A * x_in[ i,j,:] + G.Vertices[j][1].b * y_e[i,j]) .== 0
        [i in 1:G.V, j in 1:G.V; G.Adj[i][j] ≠ nothing], (G.Vertices[j][2].A * x_in[ i,j,:] + G.Vertices[j][2].b * y_e[i,j]) .== 0
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
        println(value.(y_e))
        println("The optimal path is the following:")
        println("Source s = $(s)\n")
        i = 1
        while i != 2
            j = findfirst(abs.(Vector(value.(y_e[i,:])) .-1 ) .<= 1e-4)
            print("$(Vector(value.(x_out[i,j,:]))) --> $(Vector(value.(x_in[i,j,:]))) ")
            j == 2 ? println("to reach target t") : println("at the intersection of AS $(G.Vertices[j][3]) and AS $(G.Vertices[j][4])")
            println("The cost is $(G.Adj[i][j] * sqrt(sum((Vector(value.(x_in[i,j,:])) - Vector(value.(x_out[i,j,:]))).^2)))\n")
            i = j
        end
        println("Target t = $(t)\n")
    end
end