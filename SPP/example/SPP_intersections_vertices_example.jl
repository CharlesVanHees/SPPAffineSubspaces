include("../src/AffineSubspace.jl")
include("../src/Graph.jl")
include("../src/Problem.jl")
include("../src/SPP_intersections_vertices.jl")

# Create a problem with 3 affine subspaces
A1 = AffineSubspace([0.  1. 0. 0.], [ 0.], 10.0)
A2 = AffineSubspace([3. -1. 0. 0.], [ 3.],  1.0)
A3 = AffineSubspace([3.  1. 0. 0.], [-3.],  1.0)
P = emptyProblem()
push!(P, A1)
push!(P, A2)
push!(P, A3)
setSource!(P, [-1.5, -1.5, 0., 0.])
setTarget!(P, [ 1.5, -1.5, 0., 0.])

# Create the associated graph
G = problemToGraphV(P)

# Solve the SPP problem
# The shortest path should be the following:
#    - (-3/2, -3/2) to (  0,  3  ) on AS 2   cost : 1*sqrt(3/2 ^ 2 + 9/2 ^ 2) = 3 sqrt(10) / 2
#    - ( 0  ,  3  ) to (3/2, -3/2) on AS 3   cost : 1*sqrt(3/2 ^ 2 + 9/2 ^ 2) = 3 sqrt(10) / 2
# Total cost : 3 sqrt(10)
SPPGraphV(G, P.s, P.t)