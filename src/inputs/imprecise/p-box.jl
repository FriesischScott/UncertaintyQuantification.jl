"""
	ProbabilityBox(lb::Real, up::real, dist::function, name::Symbol)

Defines an Interval, with lower a bound, an upper bound, a function to a univariate distribution from Distributions.jl and a name.
P-Box meaning:
    p_box of a Normal(μ, σ):
        lb = [μ_min, σ_min]
        up = [μ_max, σ_max]
        dist = x -> Normal(x...)

# Examples

```jldoctest
julia>  ProbabilityBox([1.75, 1.83], [1.77, 1.85], x -> Uniform(x...), :l)
ProbabilityBox(Real[1.75, 1.83], Real[1.77, 1.85], var"#7#8"(), :l)```
"""
struct ProbabilityBox{T<:UnivariateDistribution} <: ImpreciseUQInput
    lb::AbstractVector{<:Real}
    ub::AbstractVector{<:Real}
    name::Symbol

    function ProbabilityBox{T}(
        lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real}, name::Symbol
    ) where {T<:UnivariateDistribution}
        any(lb .≥ ub) && error(
            "lower bound parameters must be smaller than upper bound parameters for $name",
        )
        return new(lb, ub, name)
    end
end

function map_to_precise(
    x::Vector{<:Real}, input::ProbabilityBox{T}
) where {T<:UnivariateDistribution}
    lb = input.lb
    ub = input.ub
    name = input.name
    length(x) != length(lb) && error(
        "number of parameters $x must be equals to the number of parameter need by $name",
    )
    any(x .< lb) && error("One or more values in $x are lower than p-box's lower bound $lb")
    any(x .> ub) &&
        error("One or more values in $x are higher than p-box's upper bound $ub")

    return RandomVariable(T(x...), input.name)
end

function quantile(pbox::ProbabilityBox{T}, x) where {T<:UnivariateDistribution}

    quantiles = map(
        par -> quantile(T(par...), x),
        Iterators.product([[a, b] for (a, b) in zip(pbox.lb, pbox.ub)]...),
    )

    lb = minimum(quantiles)
    ub = maximum(quantiles)

    return Interval(lb, ub, pbox.name)
end
rand(pbox::ProbabilityBox, n::Integer=1) = quantile.(Ref(pbox), rand(n))


function sample(pbox::ProbabilityBox, n::Integer=1)
    return DataFrame(pbox.name => rand(pbox, n))
end

# Does the inverse of quantile, not cdf, which would return an interval
function reverse_quantile(pbox::ProbabilityBox{T}, x::Interval) where {T<:UnivariateDistribution}
    cdfs_lo = map(
        par -> cdf(T(par...), x.lb),
        Iterators.product([[a, b] for (a, b) in zip(pbox.lb, pbox.ub)]...),
    )

    cdfs_hi = map(
        par -> cdf(T(par...), x.ub),
        Iterators.product([[a, b] for (a, b) in zip(pbox.lb, pbox.ub)]...),
    )

    u_lo = maximum(cdfs_lo)
    u_hi = minimum(cdfs_hi)
    
    error = abs(u_hi - u_lo)
    if error >10^-6
        @warn("When inverting the quantile function for p-box $(pbox.name), the error was $(error)")
    end
    return mean([u_lo, u_hi])   # Return midpoint
end

function to_physical_space!(pbox::ProbabilityBox{T}, x::DataFrame)  where {T<:UnivariateDistribution}
    x[!, pbox.name] = _to_physical_space(pbox, x[:, pbox.name])
    return nothing
end

function _to_physical_space(d::ProbabilityBox{T}, x::Vector) where {T<:UnivariateDistribution}
    return quantile.(Ref(d), cdf.(Normal(), x))
end

function to_standard_normal_space!(pbox::ProbabilityBox{T}, x::DataFrame) where {T<:UnivariateDistribution}
    x[!, pbox.name] = _to_standard_normal_space(pbox, x[:, pbox.name])
    return nothing
end

function _to_standard_normal_space(d::ProbabilityBox{T}, x::Vector) where {T<:UnivariateDistribution}
    return quantile.(Normal(), reverse_quantile.(Ref(d), x))
end

dimensions(d::ProbabilityBox{T}) where {T<:UnivariateDistribution} = 1

length(::ProbabilityBox{T}) where {T<:UnivariateDistribution} = 1