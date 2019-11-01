struct RandomVariableSet <: AbstractInput
    members::Array{RandomVariable}
    corr::Matrix{<:Number}

    function RandomVariableSet(
        members::Array{RandomVariable},
        corr::Matrix{<:Number},
    )

        if ((length(members), length(members)) !== size(corr))
            error("wrong dimension of correlation matrix")
        end

        new(members, corr)
    end
end

# Outer constructor with default value for corr
(RandomVariableSet(
    members::Array{RandomVariable},
    corr = Matrix{Number}(I, length(members), length(members)),
) = RandomVariableSet(members, corr))

# Outer constructor for keyword passing, with default value for corr
(RandomVariableSet(
    ;
    members::Array{RandomVariable},
    corr = Matrix{Number}(I, length(members), length(members)),
) = RandomVariableSet(members, corr))

function rand(r::RandomVariableSet, n::Int64)
    # TODO: This needs to use the covariance matrix
    u = copularand(r.corr, n, length(r.members))

    samples = DataFrame()

    for (i, rv) in enumerate(r.members)
        samples[!, Symbol(rv.name)] = quantile.(rv.dist, u[:, i])
    end

    return samples
end

rand(r::RandomVariableSet) = rand(r::RandomVariableSet, 1);
