abstract type AbstractOrthogonalBasis end
struct HermiteBasis <: AbstractOrthogonalBasis
    p::Int
    d::Int
    indices::Vector{Vector{Int64}}

    HermiteBasis(p::Int, d::Int) = new(p, d, multivariate_indices(p, d))
end

struct LegendreBasis <: AbstractOrthogonalBasis
    p::Int
    d::Int
    indices::Vector{Vector{Int64}}

    LegendreBasis(p::Int, d::Int) = new(p, d, multivariate_indices(p, d))
end

struct PolynomialChaosExpansion <: UQModel
    y::Vector{Float64}
    n::Vector{Symbol}
    Ψ::AbstractOrthogonalBasis
    inputs::Vector{<:UQInput}
    output::Symbol

    function PolynomialChaosExpansion(
        data::DataFrame,
        inputs::Vector{<:UQInput},
        Ψ::AbstractOrthogonalBasis,
        output::Symbol,
    )
        random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
        random_names = names(random_inputs)

        x = data[:, random_names]
        isoprobabilistic_transform!(x, inputs, Ψ)

        A = mapreduce(row -> evaluate(collect(row), Ψ), hcat, eachrow(x))'

        y = inv(transpose(A) * A) * transpose(A) * data[:, output]

        return new(y, random_names, Ψ, random_inputs, output)
    end
end

function evaluate!(pce::PolynomialChaosExpansion, df::DataFrame)
    data = df[:, pce.n]
    isoprobabilistic_transform!(data, pce.inputs, pce.Ψ)

    out = map(row -> dot(pce.y, evaluate(collect(row), pce.Ψ)), eachrow(data))
    df[!, pce.output] = out
    return nothing
end

function sample(pce::PolynomialChaosExpansion, N::Integer)
    samps = DataFrame(randn(N, length(pce.n)), pce.n)
    out = map(row -> dot(pce.y, evaluate(collect(row), pce.Ψ)), eachrow(data))

    to_physical_space!(pce.inputs, samps)
    samps[!, pce.output] = out
    return samps
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

function He(x::Float64, d::Int)
    if d == 0
        return 1.0
    end

    if d == 1
        return x
    end

    return x * He(x, d - 1) - (d - 1) * He(x, d - 2)
end

function P(x::Float64, d::Int)
    if d == 0
        return 1.0
    end

    if d == 1
        return x
    end

    return ((2d - 1) * x * P(x, d - 1) - (d - 1) * P(x, d - 2)) / float(d)
end

function evaluate(x::Vector{Float64}, Ψ::HermiteBasis)
    return [prod(He.(x, ind)) for ind in Ψ.indices]
end

function evaluate(x::Vector{Float64}, Ψ::LegendreBasis)
    return [prod(P.(x, ind)) for ind in Ψ.indices]
end

function isoprobabilistic_transform!(
    data::DataFrame, inputs::Vector{<:UQInput}, _::HermiteBasis
)
    to_standard_normal_space!(inputs, data)
    return nothing
end

function isoprobabilistic_transform!(
    data::DataFrame, inputs::Vector{<:UQInput}, _::LegendreBasis
)
    to_standard_normal_space!(inputs, data)
    data[:, :] = quantile.(Uniform(-1, 1), cdf.(Normal(), data[:, :]))
    return nothing
end
