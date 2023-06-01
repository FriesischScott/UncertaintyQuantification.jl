function distribution_parameters(mean::Real, std::Real, _::Type{Distributions.LogNormal})
    μ = log(mean^2 / sqrt(std^2 + mean^2))
    σ = sqrt(log(std^2 / mean^2 + 1))
    return μ, σ
end
