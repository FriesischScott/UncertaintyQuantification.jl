function probabilityOfFailure(
    models::Union{Array{<:AbstractModel}, AbstractModel},
    performance::Function,
    inputs::Union{Array{<:AbstractInput}, AbstractInput},
    sim::AbstractSimulation,
)

    samples = sample(inputs, sim.n)

    # Models
    for m in models
        samples[!, Symbol(m.name)] = evaluate(m, samples)
    end

    # Probability of failure
    pf = sum(performance(samples) .< 0) / sim.n

    return pf, samples

end
