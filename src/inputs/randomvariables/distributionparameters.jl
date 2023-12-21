function distribution_parameters(mean::Real, std::Real, _::Type{Distributions.Beta})
    α = ((1 - mean) / std^2 - (1 / mean)) * mean^2
    β = α * ((1 / mean) - 1)

    return α, β
end

function distribution_parameters(mean::Real, std::Real, _::Type{Distributions.Gumbel})
    β = std * sqrt(6) / π
    μ = mean - MathConstants.γ * β
    return μ, β
end

function distribution_parameters(mean::Real, std::Real, _::Type{Distributions.LogNormal})
    μ = log(mean^2 / sqrt(std^2 + mean^2))
    σ = sqrt(log(std^2 / mean^2 + 1))
    return μ, σ
end
