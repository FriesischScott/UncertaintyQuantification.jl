struct EmpiricalDistribution <: ContinuousUnivariateDistribution
    data::Vector{<:Real}
    cdf::ECDF
    quantile::Spline1D
    pdf::Spline1D

    function EmpiricalDistribution(x::Vector{<:Real}, nbins::Integer=10)
        cdf = ecdf(x)

        f = cdf.(cdf.sorted_values)
        quantile = Spline1D(f, cdf.sorted_values)

        h = fit(Histogram, x; nbins=nbins)
        h = normalize(h; mode=:density)

        r = h.edges[1]
        pdfx = (first(r) + step(r) / 2):step(r):last(r)
        pdfy = (h.weights ./ length(x))

        pdf = Spline1D(pdfx, pdfy)
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
    return insupport(d, x) ? d.pdf(x) : zero(x)
end

function logpdf(d::EmpiricalDistribution, x::Real)
    return log(pdf(d, x))
end

function rand(rng::AbstractRNG, d::EmpiricalDistribution)
    return quantile(d, rand(rng))
end

minimum(d::EmpiricalDistribution) = minimum(d.data)
maximum(d::EmpiricalDistribution) = maximum(d.data)
mean(d::EmpiricalDistribution) = mean(d.data)
var(d::EmpiricalDistribution) = var(d.data)
