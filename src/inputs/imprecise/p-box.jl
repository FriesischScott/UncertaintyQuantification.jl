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

function sample(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution}
    x = rand()

    quantiles = map(
        par -> quantile(T(par...), x),
        Iterators.product([[a, b] for (a, b) in zip(pbox.lb, pbox.ub)]...),
    )

    lb = minimum(quantiles)
    ub = maximum(quantiles)

    return [lb, ub]
end
