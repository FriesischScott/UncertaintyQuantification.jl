struct RandomVariableSet <: RandomUQInput
    members::Array{RandomVariable}
    corr::Matrix{<:Real}

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
    corr = Matrix{Float64}(I(length(members))),
) = RandomVariableSet(members, corr))

# Outer constructor for keyword passing, with default value for corr
(RandomVariableSet(
    ;
    members::Array{RandomVariable},
    corr = Matrix{Float64}(I(length(members))),
) = RandomVariableSet(members, corr))

function sample(r::RandomVariableSet, n::Int64 = 1)
    # Sample from the gaussian copula
    ΦR = MvNormal(r.corr)
    Φ = Normal()

    u = cdf.(Φ, rand(ΦR, n))

    samples = DataFrame()

    for (i, rv) in enumerate(r.members)
        samples[!, rv.name] = quantile.(rv.dist, u[i, :])
    end

    return samples
end

function to_physical_space!(r::RandomVariableSet, x::DataFrame)
    L = cholesky(r.corr).L
    correlated_cdf = cdf.(Normal(), Matrix(x[:, names(r)]) * L')
    for (i, rv) in enumerate(r.members)
        x[!, rv.name] = quantile.(rv.dist, correlated_cdf[:, i])
    end
    return nothing
end

function to_standard_normal_space!(r::RandomVariableSet, x::DataFrame)
    for rv in r.members
        x[!, rv.name] = cdf.(rv.dist, x[:, rv.name])
    end
    L = cholesky(r.corr).L
    uncorrelated_stdnorm = quantile.(Normal(), Matrix(x[:, names(r)])) * inv(L)'
    for (i, rv) in enumerate(r.members)
        x[!, rv.name] = uncorrelated_stdnorm[:, i]
    end
    return nothing
end

names(r::RandomVariableSet) = vec(map(x -> x.name, r.members))

mean(r::RandomVariableSet) = mean(r.members)
