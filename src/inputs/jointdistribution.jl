struct JointDistribution <: RandomUQInput
    marginals::Vector{RandomVariable}
    copula::Copula

    function JointDistribution(marginals::Vector{RandomVariable}, copula::Copula)
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
        x[!, rv.name] = cdf.(rv.dist, x[:, rv.name])
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

function pdf(jd::JointDistribution, x::Vector)
    copuladensity = copuladensity(jd.copula) 

    u = Vector{eltype(x)}()
    f = 1.
    for (xi, rv) in zip(x, jd.marginals)
        push!(u, cdf(rv, xi))
        f *= pdf(rv, xi)
    end

    return copuladensity(u) * f  
end

function pdf(jd::JointDistribution, x::Matrix)
    copuladensity = copuladensity(jd.copula) 

    u = zero(x)
    f = ones(eltype(x), size(x, 1))
    for (i, rv) in enumerate(jd.marginals)
        u[:, i] = cdf.(rv.dist, x[:, i])
        f .*= pdf.(rv.dist, x[:, i])
    end

    return [copuladensity(u[i, :]) * f[i] for i in eachindex(f)] 
end
