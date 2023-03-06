function sobolindices(
    models::Vector{<:UQModel},
    inputs::Vector{<:UQInput},
    output::Symbol,
    sim::AbstractMonteCarlo,
)
    sim_double_samples = @set sim.n = 2 * sim.n

    samples = sample(inputs, sim_double_samples)
    samples = samples[shuffle(1:size(samples, 1)), :]

    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))

    evaluate!(models, samples)

    A = samples[1:(sim.n), :]
    B = samples[(sim.n + 1):end, :]

    fA = A[:, output]
    fB = B[:, output]

    fA .-= mean(fA)
    fB .-= mean(fB)

    VY = var([fA; fB])

    Si = zeros(length(random_names), 2)
    STi = zeros(length(random_names), 2)

    for (i, name) in enumerate(random_names)
        ABi = select(A, Not(name))
        ABi[:, name] = B[:, name]

        for m in models
            evaluate!(m, ABi)
        end

        ABi[:, output] .-= mean(ABi[:, output])

        first_order = x -> mean(fB .* (x .- fA)) / VY # Saltelli 2009
        total_effect = x -> (1 / (2 * sim.n)) * sum((fA .- x) .^ 2) / VY # Saltelli 2009

        # First order effects
        Si[i, 1] = first_order(ABi[:, output])
        bs = bootstrap(first_order, ABi[:, output], BasicSampling(1000))
        Si[i, 2] = stderror(bs)[1]

        # Total effects
        STi[i, 1] = total_effect(ABi[:, output])
        bs = bootstrap(total_effect, ABi[:, output], BasicSampling(1000))
        STi[i, 2] = stderror(bs)[1]
    end

    indices = DataFrame()
    indices.Variables = random_names
    indices.FirstOrder = Si[:, 1]
    indices.FirstOrderStdError = Si[:, 2]
    indices.TotalEffect = STi[:, 1]
    indices.TotalEffectStdError = STi[:, 2]

    return indices
end

function sobolindices(pce::PolynomialChaosExpansion)
    random_names = names(pce.inputs)
    indices = DataFrame()
    indices.Variables = random_names
    indices.FirstOrder =
        [
            sum(pce.y[findall(α -> α[i] != 0 && sum(α) == α[i], pce.Ψ.α)] .^ 2) for
            i in 1:length(random_names)
        ] ./ var(pce)
    indices.TotalEffect =
        [
            sum(pce.y[findall(α -> α[i] != 0, pce.Ψ.α)] .^ 2) for
            i in 1:length(random_names)
        ] ./ var(pce)

    return indices
end
