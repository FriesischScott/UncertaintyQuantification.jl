struct RandomVariable <: AbstractInput
    dist::Sampleable{Univariate}
    name::String
end

function rand(rv::RandomVariable, n::Int64 = 1)
    rand(rv.dist, n)
end
