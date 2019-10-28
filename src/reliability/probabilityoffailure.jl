function probabilityOfFailure(
    models::Array{Model},
    performance::Function,
    inputs::Array,
    sim::MonteCarlo,
)

    samples = DataFrame()

    # Parameters
    for p in filter(x -> isa(x, Parameter), inputs)
        samples[!, Symbol(p.name)] = ones(sim.n) * p.value
    end

    # RandomVariables
    for rv in filter(x -> isa(x, RandomVariable), inputs)
        samples[!, Symbol(rv.name)] = rand(rv.dist, sim.n)
    end

    # RandomVariableSets
    for rvset in filter(x -> isa(x, RandomVariableSet), inputs)
        samples = hcat(samples, rand(rvset, sim.n))
    end

    describe(samples)

    # Models
    for m in models
        samples[!, Symbol(m.name)] = evaluate(m, samples)
    end

    # Probability of failure
    sum(performance(samples) .< 0) / sim.n

end
