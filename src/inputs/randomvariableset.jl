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

function sample(r::RandomVariableSet, n::Int64 = 1)
    # Sample from the gaussian copula
    ΦR = MvNormal(r.corr)
    Φ = Normal()

    u = cdf.(Φ, rand(ΦR, n))

    samples = DataFrame()

    for (i, rv) in enumerate(r.members)
        samples[!, Symbol(rv.name)] = quantile.(rv.dist, u[i, :])
    end

    return samples
end
