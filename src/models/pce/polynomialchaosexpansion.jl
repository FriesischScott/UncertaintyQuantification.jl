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
    models::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    ls::LeastSquares,
)
    samples = sample(inputs, ls.sim)

    for m in models
        evaluate!(m, samples)
    end

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
    models::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
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

mean(pce::PolynomialChaosExpansion) = pce.y[1]
var(pce::PolynomialChaosExpansion) = sum(pce.y[2:end] .^ 2)
