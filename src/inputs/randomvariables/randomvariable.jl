"""
	RandomVariable(dist::UnivariateDistribution, name::Symbol)

Defines a random variable, with a univariate distribution from Distributions.jl and a name.

# Examples

```jldoctest
julia> RandomVariable(Normal(), :x)
RandomVariable(Normal{Float64}(μ=0.0, σ=1.0), :x)

julia> RandomVariable(Exponential(1), :x)
RandomVariable(Exponential{Float64}(θ=1.0), :x)
```
"""
struct RandomVariable{T<:Union{UnivariateDistribution,ProbabilityBox}} <: RandomUQInput
    dist::T
    name::Symbol
end

"""
	sample(rv::RandomVariable, n::Integer=1)

Generates n samples from a random variable. Returns a DataFrame.

# Examples

See also: [`RandomVariable`](@ref)
"""
function sample(rv::RandomVariable, n::Integer=1)
    return DataFrame(rv.name => rand(rv.dist, n))
end

function to_standard_normal_space!(rv::RandomVariable, x::DataFrame)
    # do nothing for standard normal rv
    if isa(rv.dist, Normal) && params(rv.dist) == (0.0, 1.0)
        return nothing
    end
    x[!, rv.name] = quantile.(Normal(), cdf.(rv.dist, x[:, rv.name]))
    return nothing
end

function to_physical_space!(rv::RandomVariable, x::DataFrame)
    # do nothing for standard normal rv
    if isa(rv.dist, Normal) && params(rv.dist) == (0.0, 1.0)
        return nothing
    end
    x[!, rv.name] = quantile.(rv.dist, cdf.(Normal(), x[:, rv.name]))
    return nothing
end

dimensions(rv::RandomVariable) = 1

logpdf(rv::RandomVariable, x::Real) = logpdf(rv.dist, x)
pdf(rv::RandomVariable, x::Real) = pdf(rv.dist, x)
cdf(rv::RandomVariable, x::Real) = cdf(rv.dist, x)
quantile(rv::RandomVariable, q::Real) = quantile(rv.dist, q)
minimum(rv::RandomVariable) = minimum(rv.dist)
maximum(rv::RandomVariable) = maximum(rv.dist)
insupport(rv::RandomVariable, x::Real) = insupport(rv.dist, x)
mean(rv::RandomVariable) = mean(rv.dist)
var(rv::RandomVariable) = var(rv.dist)
