struct RandomVariableSet
    members::Array{<:Sampleable{Univariate},2}
    names::Array{String}
    corr::Matrix{<:Number}

    function RandomVariableSet(members::Array{<:Sampleable{Univariate},2},
        names::Array,
        corr::Matrix{<:Number})


        if (length(members) !== length(names))
            error("length(members) != length(names)")
        end

        new(members, names, corr)
    end
end

# Outer constructor with default value for corr
( RandomVariableSet(members::Array{<:Sampleable{Univariate},2}, names::Array,
    corr = Matrix{Number}(I, length(members), length(members)))
    = RandomVariableSet(members,names, corr); )

# Outer constructor for keyword passing, with default value for corr
( RandomVariableSet(;members::Array{<:Sampleable{Univariate},2}, names::Array,
    corr = Matrix{Number}(I, length(members), length(members)))
    = RandomVariableSet(members,names, corr); )

function rand(r::RandomVariableSet, n::Int64)
    a = cholesky(r.corr).L
    z = rand(Normal(), n, length(r.members))
    x = cdf.(Normal(), transpose(a * transpose(z)))

    samples = DataFrame()

    for (i, (member, name)) in enumerate(zip(r.members, r.names))
        samples[Symbol(name)] = quantile.(member, x[:, i])
    end

    return samples
end

rand(r::RandomVariableSet) = rand(r::RandomVariableSet,1);
