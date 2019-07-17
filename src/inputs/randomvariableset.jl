struct RandomVariableSet
    members::Array{<:Sampleable,2}
    names::Array{String}
    corr::Matrix{Float64}

    function RandomVariableSet(members::Array{<:Sampleable,2},
        names::Array,
        corr::Matrix{Float64})


        if (length(members) !== length(names))
            error("length(members) != length(names)")
        end

        new(members, names, corr)
    end
end

# Outer constructor with default value for corr
( RandomVariableSet(members::Array{<:Sampleable,2}, names::Array,
    corr = Matrix{Float64}(I, length(members), length(members)))
    = RandomVariableSet(members,names, corr); )

# Outer constructor for keyword passing, with default value for corr
( RandomVariableSet(;members::Array{<:Sampleable,2}, names::Array,
    corr = Matrix{Float64}(I, length(members), length(members)))
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
