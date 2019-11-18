struct RandomVariable <: RandomUQInput
    dist::Sampleable{Univariate}
    name::String
end

function sample(rv::RandomVariable, n::Int64 = 1)
    DataFrame(Symbol(rv.name) => rand(rv.dist, n))
end

function to_physical_space!(rv::RandomVariable, x::DataFrame)
    x[!, Symbol(rv.name)] = quantile.(rv.dist, cdf.(Normal(), x[:, Symbol(rv.name)]))
    return nothing
end

function to_standard_normal_space!(rv::RandomVariable, x::DataFrame)
    x[!, Symbol(rv.name)] = quantile.(Normal(), cdf.(rv.dist, x[:, Symbol(rv.name)]))
    return nothing
end

mean(rv::RandomVariable) = DataFrame(Symbol(rv.name) => Distributions.mean(rv.dist))
mean(rvs::Array{RandomVariable}) = mapreduce(mean, hcat, rvs)
