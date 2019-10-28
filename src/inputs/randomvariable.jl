struct RandomVariable
    dist::Sampleable{Univariate}
    name::String
end

function rand(rv::RandomVariable, n::Int64)
    rand(rv.dist, n)
end
