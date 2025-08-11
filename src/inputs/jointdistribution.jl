struct JointDistribution{
    D<:Union{Copula,MultivariateDistribution},M<:Union{RandomVariable,Symbol}
} <: RandomUQInput
    d::D
    m::Vector{<:M}
    function JointDistribution(
        d::D, m::Vector{M}
    ) where {D<:Union{Copula,MultivariateDistribution},M<:Union{RandomVariable,Symbol}}
        if d isa Copula
            if !(eltype(m) <: RandomVariable)
                error("Must pass a marginal vector of RandomVariables.")
            end
            if dimensions(d) != length(m)
                error("Dimension mismatch between copula and marginals.")
            end
            return new{Copula,RandomVariable}(d, m)
        end
        if d isa MultivariateDistribution
            if eltype(m) != Symbol
                error("Must pass a vector of Symbols.")
            end
            if length(d) != length(m)
                error("Dimension mismatch between distribution and names.")
            end
            return new{MultivariateDistribution,Symbol}(d, m)
        end
    end
end

function sample(jd::JointDistribution{<:Copula,<:RandomVariable}, n::Integer=1)
    u = sample(jd.d, n)

    samples = DataFrame()

    for (i, rv) in enumerate(jd.m)
        samples[!, rv.name] = quantile.(rv.dist, u[:, i])
    end

    return samples
end

function to_physical_space!(jd::JointDistribution{<:Copula,<:RandomVariable}, x::DataFrame)
    correlated_cdf = to_copula_space(jd.d, Matrix{Float64}(x[:, names(jd)]))
    for (i, rv) in enumerate(jd.m)
        x[!, rv.name] = quantile.(rv.dist, correlated_cdf[:, i])
    end
    return nothing
end

function to_standard_normal_space!(
    jd::JointDistribution{<:Copula,<:RandomVariable}, x::DataFrame
)
    for rv in jd.m
        if isa(rv.dist, ProbabilityBox)
            x[!, rv.name] = reverse_quantile.(rv.dist, x[:, rv.name])
        else
            x[!, rv.name] = cdf.(rv.dist, x[:, rv.name])
        end
    end
    uncorrelated_stdnorm = to_standard_normal_space(jd.d, Matrix{Float64}(x[:, names(jd)]))
    for (i, rv) in enumerate(jd.m)
        x[!, rv.name] = uncorrelated_stdnorm[:, i]
    end
    return nothing
end

function to_standard_normal_space!(_::JointDistribution{D,M}, _::DataFrame) where {D,M}
    return error("Cannot map $D to standard normal space.")
end

function to_physical_space!(_::JointDistribution{D,M}, _::DataFrame) where {D,M}
    return error("Cannot map $D to physical space.")
end

function names(jd::JointDistribution{<:Copula,<:RandomVariable})
    return vec(map(x -> x.name, jd.m))
end

mean(jd::JointDistribution{<:Copula,<:RandomVariable}) = mean.(jd.m)

dimensions(jd::JointDistribution) = dimensions(jd.copula)

function bounds(
    jd::JointDistribution{
        <:Copula,<:RandomVariable{<:Union{UnivariateDistribution,ProbabilityBox}}
    },
)
    b = map(bounds, filter(isimprecise, jd.m))

    return vcat(getindex.(b, 1)...), vcat(getindex.(b, 2)...)
end

dimensions(jd::JointDistribution{<:Copula,<:RandomVariable}) = dimensions(jd.d)
