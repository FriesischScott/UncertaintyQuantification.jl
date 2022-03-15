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

struct LeastSquares
    sim::AbstractMonteCarlo
end

abstract type AbstractQuadrature end
struct GaussQuadrature <: AbstractQuadrature end
struct SparseQuadrature <: AbstractQuadrature end

struct PolynomialChaosExpansion <: UQModel
    y::Vector{Float64}
    Ψ::AbstractOrthogonalBasis
    output::Symbol
    inputs::Vector{<:UQInput}
end

function PolynomialChaosExpansion(
    data::DataFrame, inputs::Vector{<:UQInput}, Ψ::AbstractOrthogonalBasis, output::Symbol
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    random_names = names(random_inputs)

    x = data[:, random_names]
    isoprobabilistic_transform!(x, inputs, Ψ)

    A = mapreduce(row -> evaluate(collect(row), Ψ), hcat, eachrow(x))'

    y = inv(transpose(A) * A) * transpose(A) * data[:, output]

    return new(y, random_names, Ψ, random_inputs, output)
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    models::Vector{<:UQModel},
    Ψ::AbstractOrthogonalBasis,
    output::Symbol,
    ls::LeastSquares,
)
    samples = sample(inputs, ls.sim)

    for m in models
        evaluate!(m, samples)
    end

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    random_names = names(random_inputs)

    x = samples[:, random_names]
    isoprobabilistic_transform!(x, inputs, Ψ)

    A = mapreduce(row -> evaluate(collect(row), Ψ), hcat, eachrow(x))'
    y = inv(transpose(A) * A) * transpose(A) * samples[:, output]

    ϵ = samples[:, output] - A * y
    mse = dot(ϵ, ϵ)

    return PolynomialChaosExpansion(y, Ψ, output, inputs), samples, mse
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    models::Vector{<:UQModel},
    Ψ::AbstractOrthogonalBasis,
    output::Symbol,
    _::GaussQuadrature,
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    random_names = names(random_inputs)

    samples = DataFrame()

    if Ψ isa LegendreBasis
        _x, _w = gausslegendre(Ψ.p + 1)
    else
        _x, _w = gausshermite(Ψ.p + 1)
    end

    nodes =
        mapreduce(
            collect,
            hcat,
            Iterators.product(Iterators.repeated(_x, length(random_inputs))...),
        )'
    weights = map(prod, Iterators.product(Iterators.repeated(_w, length(random_inputs))...))

    samples = DataFrame(nodes, random_names)
    inverse_isoprobabilistic_transform!(samples, inputs, Ψ)

    for m in models
        evaluate!(m, samples)
    end

    y = zeros(length(Ψ.indices))

    for (x, w, f) in zip(eachrow(nodes), weights, samples[:, output])
        y +=
            f * w * evaluate(collect(x), Ψ) /
            (length(random_inputs) * length(random_inputs))
    end

    return PolynomialChaosExpansion(y, Ψ, output, inputs), samples
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    models::Vector{<:UQModel},
    Ψ::AbstractOrthogonalBasis,
    output::Symbol,
    _::SparseQuadrature,
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    random_names = names(random_inputs)

    samples = DataFrame()

    if Ψ isa LegendreBasis
        _x, weights = sparsegrid(length(random_inputs), Ψ.p + 1, gausslegendre)
    else
        _x, weights = sparsegrid(length(random_inputs), Ψ.p + 1, gausshermite)
    end

    nodes = reduce(hcat, _x)'

    samples = DataFrame(nodes, random_names)
    inverse_isoprobabilistic_transform!(samples, inputs, Ψ)

    for m in models
        evaluate!(m, samples)
    end

    y = zeros(length(Ψ.indices))

    for (x, w, f) in zip(eachrow(nodes), weights, samples[:, output])
        y +=
            f * w * evaluate(collect(x), Ψ) / length(random_inputs)^2
    end

    return PolynomialChaosExpansion(y, Ψ, output, inputs), samples
end

function evaluate!(pce::PolynomialChaosExpansion, df::DataFrame)
    random_inputs = filter(i -> isa(i, RandomUQInput), pce.inputs)
    random_names = names(random_inputs)

    data = df[:, random_names]
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
    normalization = [prod(1 ./ factorial.(ind)) for ind in Ψ.indices]
    return [prod(He.(x, ind)) for ind in Ψ.indices] .* normalization
end

function evaluate(x::Vector{Float64}, Ψ::LegendreBasis)
    normalization = [prod(sqrt.(2 * ind .+ 1)) for ind in Ψ.indices]
    return [prod(P.(x, ind)) for ind in Ψ.indices] .* normalization
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

function inverse_isoprobabilistic_transform!(
    data::DataFrame, inputs::Vector{<:UQInput}, _::LegendreBasis
)
    data[:, :] = quantile.(Normal(), cdf.(Uniform(-1, 1), data[:, :]))

    return to_physical_space!(inputs, data)
end

mean(pce::PolynomialChaosExpansion) = pce.y[1]
var(pce::PolynomialChaosExpansion) = sum(pce.y[2:end] .^ 2)
