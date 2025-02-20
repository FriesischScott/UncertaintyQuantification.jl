"""
    EmpiricalDistribution(x::Vector{<:Real})

    Creates an empirical distribution from the data given in `x` using kernel density estimation.
    The kernel used is Gaussian and the bandwith is obtained through the Sheather-Jones method.
    The support is inferred from the kde using numerical root finding.
"""
struct EmpiricalDistribution <: ContinuousUnivariateDistribution
    data::Vector{<:Real}
    lb::Real
    ub::Real
    h::Real
    c::Spline1D
    q::Spline1D

    function EmpiricalDistribution(data::Vector{<:Real}, n::Integer=10000)
        h = sheather_jones_bandwidth(data)

        lb = find_zero(u -> kde(h, u, data), minimum(data), Order2())

        ub = find_zero(u -> kde(h, u, data), maximum(data), Order2())

        x = collect(range(lb, ub, n))

        y = zeros(eltype(x), size(x))
        for i in eachindex(x)
            y[i] = if i == 1
                quadgk(x -> kde(h, x, data), lb, x[i])[1]
            else
                y[i - 1] + quadgk(x -> kde(h, x, data), x[i - 1], x[i])[1]
            end
        end

        clamp!(y, 0.0, 1.0)

        y[1] = 0.0
        y[end] = 1.0

        c = Spline1D(x, y; k=1, s=0.0)

        unique_idx = findlast.(isequal.(unique(y)), [y])

        q = Spline1D(y[unique_idx], x[unique_idx]; k=1, s=0.0)

        return new(data, lb, ub, h, c, q)
    end
end

function cdf(d::EmpiricalDistribution, x::Real)
    return d.c(x)
end

function quantile(d::EmpiricalDistribution, u::Real)
    return d.q(u)
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
