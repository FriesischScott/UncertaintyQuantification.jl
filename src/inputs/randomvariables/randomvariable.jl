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
struct RandomVariable{T} <: RandomUQInput where {T<:UnivariateDistribution}
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
    x[!, rv.name] = quantile.(Normal(), cdf.(rv.dist, x[:, rv.name]))
    return nothing
end

function to_physical_space!(rv::RandomVariable, x::DataFrame)
    x[!, rv.name] = quantile.(rv.dist, cdf.(Normal(), x[:, rv.name]))
    return nothing
end

function to_standard_normal_space!(rv::RandomVariable{Normal}, x::DataFrame)
    μ, σ = params(rv.dist)
    # do nothing for standard normal rv
    if μ == 0.0 && σ == 1.0
        return nothing
    else
        x[!, rv.name] = quantile.(Normal(), cdf.(rv.dist, x[:, rv.name]))
        return nothing
    end
end

function to_physical_space!(rv::RandomVariable{Normal}, x::DataFrame)
    μ, σ = params(rv.dist)
    # do nothing for standard normal rv
    if μ == 0.0 && σ == 1.0
        return nothing
    else
        x[!, rv.name] = quantile.(rv.dist, cdf.(Normal(), x[:, rv.name]))
        return nothing
    end
end

function to_standard_normal_space!(rv::RandomVariable{EmpiricalDistribution}, x::DataFrame)
    # call vectorized cdf for empirical distributions
    x[!, rv.name] = quantile.(Normal(), cdf(rv.dist, x[:, rv.name]))
    return nothing
end

function to_physical_space!(rv::RandomVariable{EmpiricalDistribution}, x::DataFrame)
    # call vectorized quantile for empirical distributions
    x[!, rv.name] = quantile(rv.dist, cdf.(Normal(), x[:, rv.name]))
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
