struct GaussianCopula <: Copula
    correlation::Matrix{<:Real}
end

const Φ = Normal()

function sample(c::GaussianCopula, n::Integer=1)
    L = cholesky(c.correlation).L

    Z = rand(Φ, n, size(c.correlation, 2))
    X = transpose(L * transpose(Z))

    return u = cdf.(Φ, X)
end

function dimensions(c::GaussianCopula)
    return size(c.correlation, 1)
end

function to_standard_normal_space(c::GaussianCopula, u)
    L = cholesky(c.correlation).L
    return quantile.(Φ, u) * inv(L)'
end

function to_copula_space(c::GaussianCopula, s)
    L = cholesky(c.correlation).L
    return cdf.(Φ, s * L')
end

function copuladensity(c::GaussianCopula)
    detC = det(c.correlation)
    M = inv(c.correlation)
    M[diagind(M)] .-= 1
    return u -> (1 / sqrt(detC)) * exp(-0.5 * quantile.(Φ, u)' * M * quantile.(Φ, u))
end