struct RandomVariableSet
    members::Array{UnivariateDistribution}
    names::Array{String}
    corr::Matrix{Number}
end

RandomVariableSet(members::Array{UnivariateDistribution}, names::Array{String}) = RandomVariableSet(members, names, Matrix{Float64}(I, 4, 4))

function rand(r::RandomVariableSet, n::Int64)
    a = cholesky(r.corr).L
    z = rand(Normal(), n, length(r.members))
    x = cdf.(Normal(), transpose(a * transpose(z)))

    samples = DataFrame()

    for (i, (member, name)) in enumerate(zip(r.members, r.names))
        samples[Symbol(name)] = quantile.(member, x[:, i])
    end

    return samples
end