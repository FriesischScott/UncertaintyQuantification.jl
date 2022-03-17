struct PolynomialChaosExpansion <: UQModel
    y::Vector{Float64}
    Ψ::PolynomialChaosBasis
    output::Symbol
    inputs::Vector{<:UQInput}
end

struct LeastSquares
    sim::AbstractMonteCarlo
end

struct GaussQuadrature end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    ls::LeastSquares,
)
    samples = sample(inputs, ls.sim)
    evaluate!(model, samples)

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    random_names = names(random_inputs)

    to_standard_normal_space!(random_inputs, samples)
    x = map_to_bases(Ψ, Matrix(samples[:, random_names]))

    A = mapreduce(row -> evaluate(Ψ, vec(row)), hcat, eachrow(x))'
    y = inv(transpose(A) * A) * transpose(A) * samples[:, output]

    ϵ = samples[:, output] - A * y
    mse = dot(ϵ, ϵ)

    to_physical_space!(random_inputs, samples)

    return PolynomialChaosExpansion(y, Ψ, output, random_inputs), samples, mse
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    _::GaussQuadrature,
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    random_names = names(random_inputs)

    samples = DataFrame()

    nodes =
        mapreduce(collect, hcat, Iterators.product(quadrature_nodes.(Ψ.p + 1, Ψ.bases)...))'
    weights = map(prod, Iterators.product(quadrature_weights.(Ψ.p + 1, Ψ.bases)...))

    samples = map_from_bases(Ψ, nodes)
    samples = DataFrame(samples, random_names)
    to_physical_space!(random_inputs, samples)

    evaluate!(model, samples)

    y = zeros(length(Ψ.α))

    for (x, w, f) in zip(eachrow(nodes), weights, samples[:, output])
        y += f * w * evaluate(Ψ, collect(x)) / length(random_inputs)^2
    end

    return PolynomialChaosExpansion(y, Ψ, output, inputs), samples
end

function evaluate!(pce::PolynomialChaosExpansion, df::DataFrame)
    data = df[:, names(pce.inputs)]
    to_standard_normal_space!(pce.inputs, data)

    data = map_to_bases(pce.Ψ, Matrix(data))

    out = map(row -> dot(pce.y, evaluate(pce.Ψ, collect(row))), eachrow(data))
    df[!, pce.output] = out
    return nothing
end

function sample(pce::PolynomialChaosExpansion, n::Integer)
    samps = map_from_bases(pce.Ψ, hcat(sample.(n, pce.Ψ.bases)...))
    out = map(row -> dot(pce.y, evaluate(pce.Ψ, collect(row))), eachrow(samps))

    samps = DataFrame(samps, names(pce.inputs))
    to_physical_space!(pce.inputs, samps)

    samps[!, pce.output] = out
    return samps
end

mean(pce::PolynomialChaosExpansion) = pce.y[1]
var(pce::PolynomialChaosExpansion) = sum(pce.y[2:end] .^ 2)
