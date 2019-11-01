function sample(inputs::Array{<:AbstractInput}, n::Int64 = 1)
    samples = DataFrame()

    # Parameters
    for p in filter(x -> isa(x, Parameter), inputs)
        samples[!, Symbol(p.name)] = ones(n) * p.value
    end

    # RandomVariables
    for rv in filter(x -> isa(x, RandomVariable), inputs)
        samples[!, Symbol(rv.name)] = rand(rv.dist, n)
    end

    # RandomVariableSets
    for rvset in filter(x -> isa(x, RandomVariableSet), inputs)
        samples = hcat(samples, rand(rvset, n))
    end

    return samples
end

sample(inputs::AbstractInput, n::Int64 = 1) = sample([inputs], n)
