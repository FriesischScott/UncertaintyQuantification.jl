function probability_of_failure(
    models::Union{Array{<:AbstractModel}, AbstractModel},
    performance::Function,
    inputs::Union{Array{<:AbstractInput}, AbstractInput},
    sim::MonteCarlo,
)

    samples = sample(inputs, sim.n)

    # Models
    for m in models
        samples = evaluate(m, samples)
    end

    # Probability of failure
    pf = sum(performance(samples) .< 0) / sim.n

    return pf, samples

end
