struct RandomVariable <: AbstractInput
    dist::Sampleable{Univariate}
    name::String
end

function sample(rv::RandomVariable, n::Int64 = 1)
    DataFrame(Symbol(rv.name) => rand(rv.dist, n))
end

mean(rv::RandomVariable) = Distributions.mean(rv.dist)
var(rv::RandomVariable) = Distributions.var(rv.dist)
