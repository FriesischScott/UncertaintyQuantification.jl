function to_standard_normal_space!(
    rv::RandomVariable{<:UnivariateDistribution}, x::DataFrame
)
    # do nothing for standard normal rv
    if isa(rv.dist, Normal) && params(rv.dist) == (0.0, 1.0)
        return nothing
    end
    x[!, rv.name] = quantile.(Normal(), cdf.(rv.dist, x[:, rv.name]))
    return nothing
end

function to_physical_space!(rv::RandomVariable{<:UnivariateDistribution}, x::DataFrame)
    # do nothing for standard normal rv
    if isa(rv.dist, Normal) && params(rv.dist) == (0.0, 1.0)
        return nothing
    end
    x[!, rv.name] = quantile.(rv.dist, cdf.(Normal(), x[:, rv.name]))
    return nothing
end
