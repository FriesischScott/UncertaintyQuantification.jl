struct EmpiricalDistribution <: ContinuousUnivariateDistribution
    data::Vector{<:Real}
    lb::Real
    ub::Real
    cdf::ECDF
    quantile::Spline1D
    h::Real

    function EmpiricalDistribution(x::Vector{<:Real}, lb::Real=-Inf, ub::Real=Inf)
        cdf = ecdf(x)

        f = cdf.(cdf.sorted_values)
        quantile = Spline1D(f, cdf.sorted_values)

        h = sheather_jones_bandwidth(x)
        return new(x, lb, ub, cdf, quantile, h)
    end
end

function cdf(d::EmpiricalDistribution, x::Real)
    return clamp(d.cdf(x), 0, 1)
end

function quantile(d::EmpiricalDistribution, x::Real)
    return d.quantile(x)
end

function pdf(d::EmpiricalDistribution, x::Real)
    return insupport(d, x) ? kde(d.h, x, d.data) : zero(x)
end

function logpdf(d::EmpiricalDistribution, x::Real)
    return log(pdf(d, x))
end

function rand(rng::AbstractRNG, d::EmpiricalDistribution)
    return quantile(d, rand(rng))
end

insupport(d::EmpiricalDistribution, x::Real) = minimum(d) <= x <= maximum(d)
minimum(d::EmpiricalDistribution) = d.lb
maximum(d::EmpiricalDistribution) = d.ub
mean(d::EmpiricalDistribution) = mean(d.data)
var(d::EmpiricalDistribution) = var(d.data)
