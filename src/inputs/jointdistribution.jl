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
    cpdf = copuladensity(jd.copula) 

    u = cdf.(jd.marginals, x)
    f = pdf.(jd.marginals, x)

    u = input_correction.(u)
    return cpdf(u) * prod(f)  
end

function pdf(jd::JointDistribution, x::DataFrame)
    cpdf = copuladensity(jd.copula) 

    cpdf_u = ones(size(x, 1))
    f = ones(size(x, 1))

    for i in axes(x, 1)
        u = cdf.(jd.marginals, Vector(x[i, names(jd.marginals)]))
        u = input_correction.(u)
        cpdf_u[i] = cpdf(u)
        
        f[i] = prod(pdf.(jd.marginals, Vector(x[i, names(jd.marginals)])))
    end
    
    return cpdf_u .* f  
end

function input_correction(u::Real)
    if u >= 1.0
        return 1.0 - eps(typeof(u))
    elseif u <= 0.0
        return 0.0 + eps(typeof(u))
    else
        return u 
    end
end
