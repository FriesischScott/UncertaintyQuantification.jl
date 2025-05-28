struct JointDistribution <: RandomUQInput
    marginals::Vector{<:RandomVariable}
    copula::Copula

    function JointDistribution(marginals::Vector{<:RandomVariable}, copula::Copula)
        length(marginals) == dimensions(copula) ||
            error("Dimension mismatch between copula and marginals")

        return new(marginals, copula)
    end
end

function sample(jd::JointDistribution, n::Integer=1)
    u = sample(jd.copula, n)

    samples = DataFrame()

    for (i, rv) in enumerate(jd.marginals)
        samples[!, rv.name] = quantile.(rv.dist, u[:, i])
    end

    return samples
end

function to_physical_space!(jd::JointDistribution, x::DataFrame)
    correlated_cdf = to_copula_space(jd.copula, Matrix{Float64}(x[:, names(jd)]))
    for (i, rv) in enumerate(jd.marginals)
        x[!, rv.name] = quantile.(rv.dist, correlated_cdf[:, i])
    end
    return nothing
end

function to_standard_normal_space!(jd::JointDistribution, x::DataFrame)
    for rv in jd.marginals
        if isa(rv.dist, ProbabilityBox)
            x[!, rv.name] = reverse_quantile.(rv.dist, x[:, rv.name])
        else
            x[!, rv.name] = cdf.(rv.dist, x[:, rv.name])
        end
    end
    uncorrelated_stdnorm = to_standard_normal_space(
        jd.copula, Matrix{Float64}(x[:, names(jd)])
    )
    for (i, rv) in enumerate(jd.marginals)
        x[!, rv.name] = uncorrelated_stdnorm[:, i]
    end
    return nothing
end

names(jd::JointDistribution) = vec(map(x -> x.name, jd.marginals))

mean(jd::JointDistribution) = mean.(jd.marginals)

dimensions(jd::JointDistribution) = dimensions(jd.copula)

function bounds(jd::JointDistribution)
    b = map(bounds, filter(isimprecise, jd.marginals))

    return vcat(getindex.(b, 1)...), vcat(getindex.(b, 2)...)
end
