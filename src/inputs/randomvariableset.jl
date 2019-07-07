struct RandomVariableSet
    members::Array{UnivariateDistribution}
    names::Array{String}
    corr::Matrix{Float64}

    function RandomVariableSet(members::Array{Distribution{Univariate, S}},
        names::Array,
        corr) where S<:ValueSupport


        if (length(members) !== length(names))
            error("length(members) != length(names)")
        end

        new(members, names, corr)
    end
end

function RandomVariableSet(members::Array{Distribution{Univariate, S}},
    names::Array) where S<:ValueSupport
    corr = Matrix{Float64}(I, length(members), length(members))

    RandomVariableSet(members, names, corr)
end

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