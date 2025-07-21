using Gurobi, JuMP

include("Graph.jl")

myOptimizer = Gurobi

function solveSPP(G::DirectedGraph, s::T, t::T, verbose::Bool=true) where {T}
    @assert T <: eltype(G.Vertices)
    @assert s in G.Vertices && t in G.Vertices

    model = Model(myOptimizer.Optimizer)
    set_silent(model)
    
    @variable(model, y_e[1:G.V, 1:G.V], Bin)
    for i in 1:G.V, j in 1:G.V
        if isnothing(G.Adj[i][j]) fix(y_e[i,j], 0) end # no edge between i and j
    end

    @objective(model, Min, sum(G.Adj[i][j] ≠ nothing ? y_e[i,j] * G.Adj[i][j] : 0 for j in 1:G.V, i in 1:G.V))

    # Flow conservation constraint
    @constraint(model, [i in 1:G.V], (sum(y_e[:,i]) + (s == G.Vertices[i])) == (sum(y_e[i,:]) + (t == G.Vertices[i])))
    # Degree constraint
    @constraint(model, [i in 1:G.V], (sum(y_e[:,i])) + (s == G.Vertices[i]) <= 1)

    optimize!(model); println()
    if is_solved_and_feasible(model)
        println("A solution of optimal cost $(objective_value(model)) has been found!")
    else println("There is no path between your source and target in the given graph."); return
    end

    if verbose
        println("The optimal path is the following:")
        println("Source s = $(s)")
        i = findfirst(G.Vertices .== s)
        while !(i in findall(G.Vertices .== t))
            j = findfirst(Vector(value.(y_e[i,:])) .== 1)
            println("$(G.Vertices[i]) --> $(G.Vertices[j])\tCost: $(G.Adj[i][j])")
            i = j
        end
        println("Target t = $(t)")
    end
end