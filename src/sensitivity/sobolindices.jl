function sobolindices(
    models::Array{<:UQModel},
    inputs::Array{<:UQInput},
    output::Symbol,
    sim::AbstractMonteCarloSampling
)

    sim_double_samples = @set sim.n = 2 * sim.n

    samples = sample(inputs, sim_double_samples)

    random_names = filter(i -> isa(i, RandomUQInput), inputs) |> names

    for m ∈ models
        evaluate!(m, samples)
    end

    A = samples[1:sim.n, :]
    B = samples[sim.n + 1:end, :]

    fA = A[:, output]
    fB = B[:, output]

    fA .-= mean(fA)
    fB .-= mean(fB)

    VY = var([fA; fB])

    Si = zeros(length(random_names), 2)
    STi = zeros(length(random_names), 2)

    for (i, name) ∈ enumerate(random_names)
        ABi = select(A, Not(name))
        ABi[:, name] = B[:, name]

        for m ∈ models
            evaluate!(m, ABi)
        end

        ABi[:, output] .-= mean(ABi[:, output])

        first_order = x -> mean(fB .* (x .- fA )) / VY # Saltelli 2009
        total_effect = x ->  (1 / (2 * sim.n)) * sum((fA .- x).^2) / VY # Saltelli 2009
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
