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

    VY = var([fA; fB])
    Si = zeros(length(random_names))
    STi = zeros(length(random_names))

    for (i, name) ∈ enumerate(random_names)
        ABi = select(A, Not(name))
        ABi[:, name] = B[:, name]

        for m ∈ models
            evaluate!(m, ABi)
        end

        # First order effects
        Si[i] = mean(fB .* (ABi[:, output] .- fA )) / VY # Saltelli 2009
        # Total effects
        STi[i] = (1 / (2 * sim.n)) * sum((fA .- ABi[:, output]).^2) / VY # Saltelli 2009
    end

    return (;zip(random_names, Si)...), (;zip(random_names, STi)...)
end
