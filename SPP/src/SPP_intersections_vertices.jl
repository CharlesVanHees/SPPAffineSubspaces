using Gurobi, JuMP

include("AffineSubspace.jl")
include("Graph.jl")
include("Problem.jl")

myOptimizer = Gurobi

function problemToGraphV(P::Problem)
    
end