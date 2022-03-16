abstract type AbstractOrthogonalBasis end

struct PolynomialChaosBasis
    bases::Vector{<:AbstractOrthogonalBasis}
    p::Int
    d::Int
    α::Vector{Vector{Int64}}

    function PolynomialChaosBasis(bases::Vector{<:AbstractOrthogonalBasis}, p::Int)
        d = length(bases)
        return new(bases, p, d, multivariate_indices(p, d))
    end
end

function evaluate(Ψ::PolynomialChaosBasis, x::AbstractVector{Float64})
    return float.(Symbolics.value.([prod(evaluate.(Ψ.bases, x, α)) for α in Ψ.α]))
end

struct LegendreBasis <: AbstractOrthogonalBasis
    p::Int
    P::Vector{Num}

    LegendreBasis(p::Int, normalize::Bool=true) = new(p, legendre(p, normalize))
end

function evaluate(Ψ::LegendreBasis, x::Real, d::Int)
    return Ψ.P[d + 1](x)
end

function legendre(p::Int, normalize::Bool=true)
    @variables x

    P = Vector(undef, p + 1)
    P[1] = 1
    P[2] = x

    for n in 2:(p)
        P[n + 1] = ((2n - 1) * x * P[n] - (n - 1) * P[n - 1]) / n
    end

    if normalize
        P[2] = P[2] * sqrt(3)

        for n in 2:(p)
            P[n + 1] = P[n + 1] * sqrt(2n + 1)
        end
    end

    return build_function.(P, x, expression=false)
end

function multivariate_indices(p::Int, d::Int)
    No = Int64(factorial(p + d) / factorial(p) / factorial(d))

    idx = vcat(zeros(Int64, 1, d), Matrix(I, d, d), zeros(Int64, No - d - 1, d))

    pᵢ = ones(Int64, d, No)

    for k in 2:No
        for i in 1:d
            pᵢ[i, k] = sum(pᵢ[i:d, k - 1])
        end
    end

    P = d + 1
    for k in 2:p
        L = P
        for j in 1:d, m in (L - pᵢ[j, k] + 1):L
            P += 1
            idx[P, :] = idx[m, :]
            idx[P, j] = idx[P, j] + 1
        end
    end

    return map(collect, eachrow(idx))
end

function map_to_base(_::LegendreBasis, x::AbstractVector)
    return quantile.(Uniform(-1, 1), cdf.(Normal(), x))
end

function map_to_bases(Ψ::PolynomialChaosBasis, x::AbstractMatrix)
    return hcat([map_to_base(Ψ.bases[i], x[:, i]) for i in 1:length(Ψ.bases)]...)
end
