struct EmpiricalDistribution <: ContinuousUnivariateDistribution
    data::Vector{<:Real}
    cdf::ECDF
    quantile::Spline1D
    pdf::InterpKDE

    function EmpiricalDistribution(x::Vector{<:Real}, kernel=Normal)
        cdf = ecdf(x)

        f = cdf.(cdf.sorted_values)
        quantile = Spline1D(f, cdf.sorted_values)

        pdf = InterpKDE(kde_lscv(x; kernel=kernel))
        return new(x, cdf, quantile, pdf)
    end
end

function cdf(d::EmpiricalDistribution, x::Real)
    return clamp(d.cdf(x), 0, 1)
end

function quantile(d::EmpiricalDistribution, x::Real)
    return d.quantile(x)
end

function pdf(d::EmpiricalDistribution, x::Real)
    return insupport(d, x) ? abs(pdf(d.pdf, x)) : zero(x)
end

function logpdf(d::EmpiricalDistribution, x::Real)
    return log(pdf(d, x))
end

function rand(rng::AbstractRNG, d::EmpiricalDistribution)
    return quantile(d, rand(rng))
end

insupport(d::EmpiricalDistribution, x::Real) = minimum(d) <= x <= maximum(d)
minimum(d::EmpiricalDistribution) = minimum(d.data)
maximum(d::EmpiricalDistribution) = maximum(d.data)
mean(d::EmpiricalDistribution) = mean(d.data)
var(d::EmpiricalDistribution) = var(d.data)
