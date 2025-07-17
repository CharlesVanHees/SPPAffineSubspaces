include("./AffineSubspace.jl")

mutable struct Problem
    M::Int64 # Number of affine subspaces

    affSubspaces::Vector{AffineSubspace} # List of the affine subspaces
    intersections::Vector{Vector{UInt64}} # For each affine subspace, list of the index of the affine subspaces it intersects.

    s::Vector{Real} # source
    t::Vector{Real} # target
end

emptyProblem() = Problem(0, AffineSubspace[],  Vector{UInt64}[], Real[], Real[])

setSource!(P::Problem, s::Vector{T}) where T <: Real  = P.s = s;
setTarget!(P::Problem, t::Vector{T}) where T <: Real = P.t = t;

function Base.:push!(P::Problem, AS::AffineSubspace)
    P.M += 1
    push!(P.intersections, UInt64[])
    for (i, as) in pairs(P.affSubspaces)
        if intersect(as, AS)
            push!(P.intersections[i], P.M)
            push!(P.intersections[P.M], i)
        end
    end
    push!(P.affSubspaces, AS)
end

function printProblem(P::Problem)
    println("Number of affine subspaces: $(P.M)")
    println("Those are:")
    for AS in P.affSubspaces println(AS) end
    println("\nIntersection lists:")
    for i in 1:P.M
        println("Affine subspace nb $(i): $(Vector{Int64}(P.intersections[i]))")
    end
end

function exampleProblem()
    A1 = AffineSubspace([0 1], [0], 10.0)
    A2 = AffineSubspace([3 -1], [3], 1.0)
    A3 = AffineSubspace([3 1], [-3], 1.0)
    A4 = AffineSubspace([0 1], [-1], 4.0)
    P = emptyProblem()
    push!(P, A1)
    push!(P, A2)
    push!(P, A3)
    push!(P, A4)

    setSource!(P, [-1.5, -1.5])
    setTarget!(P, [1.5, -1.5])

    return P
end

# printProblem(example())