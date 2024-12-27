struct EmpiricalDistribution <: ContinuousUnivariateDistribution
    data::Vector{<:Real}
    lb::Real
    ub::Real
    h::Real

    function EmpiricalDistribution(x::Vector{<:Real}, lb::Real=-Inf, ub::Real=Inf)
        h = sheather_jones_bandwidth(x)

        lb = find_zero(u -> kde(h, u, x), minimum(x), Order2())

        ub = find_zero(u -> kde(h, u, x), maximum(x), Order2())

        return new(x, lb, ub, h)
    end
end

function cdf(d::EmpiricalDistribution, x::Real)
    return quadgk(x -> pdf(d, x), d.lb, x)[1]
end

function quantile(d::EmpiricalDistribution, u::Real)
    return find_zero(x -> cdf(d, x) - u, (d.lb, d.ub))
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

insupport(d::EmpiricalDistribution, x::Real) = d.lb <= x <= d.ub
minimum(d::EmpiricalDistribution) = d.lb
maximum(d::EmpiricalDistribution) = d.ub
mean(d::EmpiricalDistribution) = mean(d.data)
var(d::EmpiricalDistribution) = var(d.data)
