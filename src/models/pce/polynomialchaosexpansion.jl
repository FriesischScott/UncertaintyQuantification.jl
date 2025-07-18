struct PolynomialChaosExpansion <: UQModel
    y::Vector{Float64}
    Ψ::PolynomialChaosBasis
    output::Symbol
    inputs::Vector{<:UQInput}
end

struct LeastSquares
    sim::AbstractMonteCarlo
end

"""
    WeightedApproximateFetekePoints(sim::AbstractMonteCarlo; fadd=10, fmult=2)

Struct for performing weighted approximate Feteke points (wafp) subsampling of a Monte-Carlo sampler for use in generating a 
`PolynomialChaosExpansion`. Given a `PolynomialChaosBasis` of dimension `N`, and a Monte-Carlo sampler with `M` samples, generates
a subsample of size `max(N,min(N+fadd,N+fmult,M))` biased towards maximizing the determinant of the Gramian typically requiring less 
than `M` model evaluations. Follows procedure described in [burkEfficientSampling](@cite). 
"""
struct WeightedApproximateFetekePoints
    sim::AbstractMonteCarlo
    fadd::Integer
    fmult::Integer
    function WeightedApproximateFetekePoints(sim::AbstractMonteCarlo; fadd=10, fmult=2)
        return new(sim, fadd, fmult)
    end
end

struct GaussQuadrature end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::Vector{<:UQModel},
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
    y = A \ samples[:, output]

    ϵ = samples[:, output] - A * y
    mse = mean(ϵ .^ 2)

    to_physical_space!(random_inputs, samples)

    return PolynomialChaosExpansion(y, Ψ, output, random_inputs), samples, mse
end

function polynomialchaos(
    inputs::UQInput,
    model::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    ls::LeastSquares,
)
    return polynomialchaos([inputs], model, Ψ, output, ls)
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    ls::LeastSquares,
)
    return polynomialchaos(inputs, [model], Ψ, output, ls)
end

function polynomialchaos(
    inputs::UQInput,
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    ls::LeastSquares,
)
    return polynomialchaos([inputs], [model], Ψ, output, ls)
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    wafp::WeightedApproximateFetekePoints
)
    samples = sample(inputs, wafp.sim)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    random_names = names(random_inputs)
    to_standard_normal_space!(random_inputs, samples)
    x = map_to_bases(Ψ, Matrix(samples[:, random_names]))

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)


    Np = length(Ψ.α)
    n = wafp.sim.n
    rest = min((Np-1) * wafp.fmult, wafp.fadd, n - Np)
    rest = max(rest, 0)
    
    A = Matrix{Float64}(undef, n, Np)
    for i in 1:n
        A[i,:] .= evaluate(Ψ, x[i,:])
    end
    w = norm.(eachrow(A)) .^ (-2.0)
    B = A' .* reshape(w .^ (1/2), 1, :)
    _, _, p = qr(B, ColumnNorm())
    pout = zeros(Int, Np + rest)
    pout[1:Np] .= p[1:Np]
    Ginv = inv(B[:,p] * B[:,p]')
    for i in 1:rest
        val, j = findmax(j -> B[:,p[j]]' * Ginv * B[:,p[j]], Np+i:n)
        pout[Np+i] = p[j]
        Ginv .-= ((Ginv * B[:,p[j]]) * (B[:,p[j]]' * Ginv)) ./ (1 + val)
    end
    
    samples = samples[pout,:]
    to_physical_space!(random_inputs, samples)
    w = w[pout]
    evaluate!(model, samples)

    A = A[pout,:]
    W = Diagonal(w)
    y = (A' * W * A) \ (A' * W * samples[:, output])

    ϵ = samples[:, output] - A * y
    mse = mean(ϵ .^ 2)

    return PolynomialChaosExpansion(y, Ψ, output, random_inputs), samples, mse
end

function polynomialchaos(
    inputs::UQInput,
    model::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    wafp::WeightedApproximateFetekePoints
)
    return polynomialchaos([inputs], model, Ψ, output, wafp)
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    wafp::WeightedApproximateFetekePoints
)
    return polynomialchaos(inputs, [model], Ψ, output, wafp)
end

function polynomialchaos(
    inputs::UQInput,
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    wafp::WeightedApproximateFetekePoints,
)
    return polynomialchaos([inputs], [model], Ψ, output, wafp)
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    _::GaussQuadrature,
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)
    random_names = names(random_inputs)

    nodes = mapreduce(
        n -> [n...]', vcat, Iterators.product(quadrature_nodes.(Ψ.p + 1, Ψ.bases)...)
    )
    weights = map(prod, Iterators.product(quadrature_weights.(Ψ.p + 1, Ψ.bases)...))

    samples = DataFrame(map_from_bases(Ψ, nodes), random_names)
    to_physical_space!(random_inputs, samples)

    if !isempty(deterministic_inputs)
        DataFrames.hcat!(samples, sample(deterministic_inputs, size(nodes, 1)))
    end

    evaluate!(model, samples)

    y = mapreduce(
        (x, w, f) -> f * w * evaluate(Ψ, collect(x)),
        +,
        eachrow(nodes),
        weights,
        samples[:, output],
    )

    return PolynomialChaosExpansion(y, Ψ, output, random_inputs), samples
end

function polynomialchaos(
    inputs::UQInput,
    model::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    gq::GaussQuadrature,
)
    return polynomialchaos([inputs], model, Ψ, output, gq)
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    gq::GaussQuadrature,
)
    return polynomialchaos(inputs, [model], Ψ, output, gq)
end

function polynomialchaos(
    inputs::UQInput,
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    gq::GaussQuadrature,
)
    return polynomialchaos([inputs], [model], Ψ, output, gq)
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
    samps = hcat(sample.(n, pce.Ψ.bases)...)
    out = map(row -> dot(pce.y, evaluate(pce.Ψ, collect(row))), eachrow(samps))

    samps = DataFrame(map_from_bases(pce.Ψ, samps), names(pce.inputs))
    to_physical_space!(pce.inputs, samps)

    samps[!, pce.output] = out
    return samps
end

mean(pce::PolynomialChaosExpansion) = pce.y[1]
var(pce::PolynomialChaosExpansion) = sum(pce.y[2:end] .^ 2)
