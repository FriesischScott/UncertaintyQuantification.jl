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
struct ProbabilityBox <: ImpreciseUQInput
    lb::AbstractVector{Real}
    ub::AbstractVector{Real}
    dist::Function
    name::Symbol
end


function map_to_precise(x::Vector{<:Real}, input::ProbabilityBox)
    length(x) != length(input.lb) && error("number of parameters $x must be equals to the number of parameter need by $input.name")
    any(x .< input.lb) && error("One or more values in $x are lower than p-box's lower bound $input.lb")
    any(x .> input.ub) && error("One or more values in $x are higher than p-box's upper bound $input.lb")

    return RandomVariable(input.dist(x), input.name)
end