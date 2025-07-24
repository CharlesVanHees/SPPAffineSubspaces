using Gurobi, JuMP

struct AffineSubspace
    A::Matrix # Matrix of dimension d x n
    b::Vector # Vector of size d
    β::Number    # Cost of the affine subspace
end

Base.:(==)(AS1::AffineSubspace, AS2::AffineSubspace) = AS1.A == AS2.A && AS1.b == AS2.b && AS1.β == AS2.β

function intersect(AS1::AffineSubspace, AS2::AffineSubspace; Optimizer::Module = Gurobi)
    d1 = size(AS1.A, 1); d2 = size(AS2.A, 1)
    @assert (n = size(AS1.A, 2)) == size(AS2.A, 2)

    redirect_stdout(devnull) do
        model = Model(Optimizer.Optimizer)
        set_silent(model)
        @variable(model, x[1:n])
        @constraint(model, AS1.A * x .== AS1.b)
        @constraint(model, AS2.A * x .== AS2.b)
        optimize!(model)
        return is_solved_and_feasible(model)
    end
end

contains(AS::AffineSubspace, x::Vector) = sqrt(sum((AS.A * x + AS.b).^2)) <= 1e-5