function probabilityOfFailure(
    models::Array{Model},
    performance::Function,
    inputs::Array,
    sim::AbstractSimulation,
)

    samples = sample(inputs, sim.n)

    # Models
    for m in models
        samples[!, Symbol(m.name)] = evaluate(m, samples)
    end

    # Probability of failure
    sum(performance(samples) .< 0) / sim.n

end
