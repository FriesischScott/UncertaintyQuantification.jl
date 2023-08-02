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
struct RandomVariable <: RandomUQInput
    dist::UnivariateDistribution
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
    x[!, rv.name] = _to_physical_space(rv.dist, x[:, rv.name])
    return nothing
end

function _to_physical_space(d::UnivariateDistribution, x::Vector)
    return quantile.(d, cdf.(Normal(), x))
end

function _to_physical_space(d::Normal, x::Vector)
    if d.μ == 0.0 && d.σ == 1.0
        return x
    else
        return x .* d.σ .+ d.μ
    end
end

function _to_physical_space(d::LogNormal, x::Vector)
    return exp.(x .* d.σ .+ d.μ)
end

function _to_physical_space(d::Uniform, x::Vector)
    return cdf.(Normal(), x) .* (d.b - d.a) .+ d.a
end

function to_standard_normal_space!(rv::RandomVariable, x::DataFrame)
    x[!, rv.name] = _to_standard_normal_space(rv.dist, x[:, rv.name])
    return nothing
end

function _to_standard_normal_space(d::UnivariateDistribution, x::Vector)
    return quantile.(Normal(), cdf.(d, x))
end

function _to_standard_normal_space(d::Normal, x::Vector)
    if d.μ == 0.0 && d.σ == 1.0
        return x
    else
        return (x .- d.μ) ./ d.σ
    end
end

function _to_standard_normal_space(d::LogNormal, x::Vector)
    return (log.(x) .- d.μ) ./ d.σ
end

function _to_standard_normal_space(d::Uniform, x::Vector)
    return quantile.(Normal(), (x .- d.a) ./ (d.b - d.a))
end

dimensions(rv::RandomVariable) = 1

logpdf(rv::RandomVariable, x::Real) = logpdf(rv.dist, x)
pdf(rv::RandomVariable, x::Real) = pdf(rv.dist, x)
pdf(rv::RandomVariable, x::DataFrame) = pdf.(rv.dist, x[!, rv.name])
cdf(rv::RandomVariable, x::Real) = cdf(rv.dist, x)
quantile(rv::RandomVariable, q::Real) = quantile(rv.dist, q)
minimum(rv::RandomVariable) = minimum(rv.dist)
maximum(rv::RandomVariable) = maximum(rv.dist)
insupport(rv::RandomVariable, x::Real) = insupport(rv.dist, x)
mean(rv::RandomVariable) = mean(rv.dist)
var(rv::RandomVariable) = var(rv.dist)

function pdf(inputs::Vector{<:UQInput}, x::DataFrame)
    return mapreduce(i -> pdf(i, x), hcat, inputs)
end