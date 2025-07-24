include("./AffineSubspace.jl")

mutable struct Problem
    M::Int64 # Number of affine subspaces

    affSubspaces::Vector{AffineSubspace} # List of the affine subspaces
    intersections::Vector{Vector{UInt64}} # For each affine subspace, list of the index of the affine subspaces it intersects.

    s::Vector # source
    containSource::Vector{Bool} # vector with the same shape as affSubspaces, with true for the affine subspaces that contain the source
    t::Vector # target
    containTarget::Vector{Bool} # vector with the same shape as affSubspaces, with true for the affine subspaces that contain the target
end

emptyProblem() = Problem(0, AffineSubspace[],  Vector{UInt64}[], Real[], Bool[], Real[], Bool[])

function setSource!(P::Problem, s::Vector)
    P.s = s
    for i in eachindex(P.affSubspaces) P.containSource[i] = contains(P.affSubspaces[i], s) end
end
function setTarget!(P::Problem, t::Vector)
    P.t = t
    for i in eachindex(P.affSubspaces) P.containTarget[i] = contains(P.affSubspaces[i], t) end
end

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

    push!(P.containSource, false)
    if (!isempty(P.s) && contains(AS, P.s)) P.containSource[end] = true end
    push!(P.containTarget, false)
    if (!isempty(P.t) && contains(AS, P.t)) P.containTarget[end] = true end
end

function printProblem(P::Problem)
    println("Number of affine subspaces: $(P.M)")
    println("Those are:")
    for AS in P.affSubspaces println(AS) end
    println("\nIntersection lists:")
    for i in 1:P.M
        println("Affine subspace nb $(i): $(Vector{Int64}(P.intersections[i]))")
    end
    println("\nThe source is $(P.s).")
    print("It is contained in the affine subspaces number ")
    for i in eachindex(P.containSource)
        if P.containSource[i] print("$(i) ") end
    end
    println(".")
    println("The target is $(P.t)")
    print("It is contained in the affine subspaces number ")
    for i in eachindex(P.containTarget)
        if P.containTarget[i] print("$(i) ") end
    end
    println(".\n")
end