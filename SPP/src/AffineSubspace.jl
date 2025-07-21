using LazySets

struct AffineSubspace{T <: Real}
    A::Matrix{T} # Matrix of dimension d x n
    b::Vector{T} # Vector of size d
    β::Number    # Cost of the affine subspace
end

Base.:(==)(AS1::AffineSubspace, AS2::AffineSubspace) = AS1.A == AS2.A && AS1.b == AS2.b && AS1.β == AS2.β

function intersect(AS1::AffineSubspace, AS2::AffineSubspace)
    d1 = size(AS1.A, 1); d2 = size(AS2.A, 1)
    @assert (n = size(AS1.A, 2)) == size(AS2.A, 2)

    H = Hyperplane(Vector{Float64}(AS1.A[1,:]), -Float64(AS1.b[1]))
    for i in 2:d1 H = H ∩ Hyperplane(Vector{Float64}(AS1.A[i,:]), -Float64(AS1.b[i])) end
    for i in 1:d2 H = H ∩ Hyperplane(Vector{Float64}(AS2.A[i,:]), -Float64(AS2.b[i])) end

    return !isempty(H)
end

contains(AS::AffineSubspace, x::Vector{T}) where {T <: Real} = norm(AS.A * x + AS.b) <= 1e-8