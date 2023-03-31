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
    return [prod(evaluate.(Ψ.bases, x, α)) for α in Ψ.α]
end

struct LegendreBasis <: AbstractOrthogonalBasis
    p::Int
    P::Vector{Num}

    LegendreBasis(p::Int, normalize::Bool=true) = new(p, legendre(p, normalize))
end

struct HermiteBasis <: AbstractOrthogonalBasis
    normalize::Bool
end

HermiteBasis() = HermiteBasis(true)

function evaluate(Ψ::AbstractOrthogonalBasis, x::AbstractVector{<:Real}, d::Int)
    return map(xᵢ -> evaluate(Ψ, xᵢ, d), x)
end

function evaluate(Ψ::LegendreBasis, x::Real, d::Int)
    return Ψ.P[d + 1](x).val
end

function evaluate(Ψ::HermiteBasis, x::Real, d::Int)
    val = He(x, d)
    return Ψ.normalize ? val / sqrt(factorial(d)) : val
end

function legendre(p::Int, normalize::Bool=true)
    @variables x

    P = Vector{Num}(undef, p + 1)
    P[1] = 1.0
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

    return eval.(build_function.(P, x))
end

function He(x::Real, n::Integer)
    if n == 0
        return 1.0
    elseif n == 1
        return x
    else
        He⁻ = 1.0
        He = x

        He⁺ = 0.0
        for i in 2:n
            He⁺ = x * He - (i - 1) * He⁻
            He⁻ = He
            He = He⁺
        end
        return He⁺
    end
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

function map_to_base(_::HermiteBasis, x::AbstractVector)
    return x
end

function map_to_bases(Ψ::PolynomialChaosBasis, x::AbstractMatrix)
    return mapreduce((b, xᵢ) -> map_to_base(b, xᵢ), hcat, Ψ.bases, eachcol(x))
end

function map_from_base(_::LegendreBasis, x::AbstractVector)
    return quantile.(Normal(), cdf.(Uniform(-1, 1), x))
end

function map_from_base(_::HermiteBasis, x::AbstractVector)
    return x
end

function map_from_bases(Ψ::PolynomialChaosBasis, x::AbstractMatrix)
    return mapreduce((b, xᵢ) -> map_from_base(b, xᵢ), hcat, Ψ.bases, eachcol(x))
end

function quadrature_nodes(n::Int, _::LegendreBasis)
    x, _ = gausslegendre(n)
    return x
end

function quadrature_weights(n::Int, _::LegendreBasis)
    _, w = gausslegendre(n)
    return w ./ 2
end

function quadrature_nodes(n::Int, _::HermiteBasis)
    x, _ = gausshermite(n)
    return sqrt(2) * x
end

function quadrature_weights(n::Int, _::HermiteBasis)
    _, w = gausshermite(n)
    return w / sqrt(π)
end

function sample(n::Int, _::LegendreBasis)
    return rand(Uniform(-1, 1), n)
end

function sample(n::Int, _::HermiteBasis)
    return rand(Normal(), n)
end
