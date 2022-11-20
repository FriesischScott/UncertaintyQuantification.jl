"""
	RandomVariable(dist::Sampleable{Univariate}, name::Symbol)

Defines a random variable, with a univariate distribution from Distributions.jl and a name.

# Examples

```jldoctest
julia> RandomVariable(Normal(), :x)
RandomVariable(Normal{Float64}(μ=0.0, σ=1.0), :x)

julia> RandomVariable(Exponential(1), :x)
RandomVariable(Exponential{Float64}(θ=1.0), :x)
```
"""
struct RandomVariable <: RandomUQInput
    dist::Sampleable{Univariate}
    name::Symbol
end

"""
	sample(rv::RandomVariable, n::Integer)

Generates n samples from a random variable. Returns a DataFrame.

# Examples



See also: [`RandomVariable`](@ref)
"""
function sample(rv::RandomVariable, n::Integer=1)
    return DataFrame(rv.name => rand(rv.dist, n))
end

function to_physical_space!(rv::RandomVariable, x::DataFrame)
    x[!, rv.name] = quantile.(rv.dist, cdf.(Normal(), x[:, rv.name]))
    return nothing
end

function to_standard_normal_space!(rv::RandomVariable, x::DataFrame)
    x[!, rv.name] = quantile.(Normal(), cdf.(rv.dist, x[:, rv.name]))
    return nothing
end

mean(rv::RandomVariable) = DataFrame(rv.name => Distributions.mean(rv.dist))
mean(rvs::Array{RandomVariable}) = mapreduce(mean, hcat, rvs)

dimensions(rv::RandomVariable) = 1
