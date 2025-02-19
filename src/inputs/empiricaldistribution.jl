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

    function EmpiricalDistribution(x::Vector{<:Real})
        h = sheather_jones_bandwidth(x)

        lb = find_zero(u -> kde(h, u, x), minimum(x), Order2())

        ub = find_zero(u -> kde(h, u, x), maximum(x), Order2())

        return new(x, lb, ub, h)
    end
end

function cdf(d::EmpiricalDistribution, x::Real)
    return quadgk(x -> pdf(d, x), d.lb, x)[1]
end

# vectorized cdf function exploiting the monotonicity
function cdf(d::EmpiricalDistribution, x::AbstractVector{<:Real})
    u = zeros(eltype(x), size(x))
    idx = sortperm(x)
    for (i, xᵢ) in enumerate(idx)
        u[xᵢ] = if i == 1
            quadgk(x -> pdf(d, x), d.lb, x[xᵢ])[1]
        else
            u[idx[i - 1]] + quadgk(x -> pdf(d, x), x[idx[i - 1]], x[xᵢ])[1]
        end
    end
    return u
end

function quantile(d::EmpiricalDistribution, u::Real)
    return find_zero(x -> cdf(d, x) - u, (d.lb, d.ub), Roots.A42())
end

# vectorized quantile function exploiting the monotonicity
function quantile(d::EmpiricalDistribution, u::AbstractVector{<:Real})
    x = zeros(eltype(u), size(u))
    idx = sortperm(u)
    for (i, uᵢ) in enumerate(idx)
        x[uᵢ] = find_zero(
            x -> cdf(d, x) - u[uᵢ], (i == 1 ? d.lb : x[idx[i - 1]], d.ub), Roots.A42()
        )
    end
    return x
end

function pdf(d::EmpiricalDistribution, x::Real)
    return insupport(d, x) ? kde(d.h, x, d.data) : zero(x)
end

function logpdf(d::EmpiricalDistribution, x::Real)
    return log(pdf(d, x))
end

function Distributions._rand!(
    rng::AbstractRNG, d::EmpiricalDistribution, A::AbstractArray{<:Real}
)
    A[:] = quantile(d, rand(rng, length(A)))
    return A
end

insupport(d::EmpiricalDistribution, x::Real) = d.lb <= x <= d.ub
minimum(d::EmpiricalDistribution) = d.lb
maximum(d::EmpiricalDistribution) = d.ub
mean(d::EmpiricalDistribution) = mean(d.data)
var(d::EmpiricalDistribution) = var(d.data)
