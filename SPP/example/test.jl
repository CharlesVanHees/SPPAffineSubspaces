include("../src/AffineSubspace.jl")
include("../src/Graph.jl")
include("../src/Problem.jl")
include("../src/SPP_intersections_edges.jl")
include("../src/SPP_intersections_vertices.jl")

n = 4

P = emptyProblem()

for i in 1:8
    m = rand(1:n-1); push!(P, AffineSubspace(rand([-1,1], (m,n)) .* rand((m,n)), rand([-1,1], m) .* rand(m), rand(1:5)))
end

push!(P, AffineSubspace(zeros((0,n)), zeros(0), 100))
push!(P, AffineSubspace(zeros((0,n)), zeros(0), 100))

setSource!(P, rand([-1,1], n) .* rand(n))
setTarget!(P, rand([-1,1], n) .* rand(n))

printProblem(P)

G_E = problemToGraphE(P)
G_V = problemToGraphV(P)

SPPGraphE(G_E, P.s, P.t, R=10, verbose=true)
SPPGraphV(G_V, P.s, P.t, R=10, verbose=true)