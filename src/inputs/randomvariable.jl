struct RandomVariable <: RandomUQInput
    dist::Sampleable{Univariate}
    name::Symbol
end

function sample(rv::RandomVariable, n::Int64 = 1)
    DataFrame(rv.name => rand(rv.dist, n))
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
