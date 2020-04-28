struct GaussianCopula <: Copula
    correlation::Matrix{<:Real}
end

const Φ = Normal()

function sample(c::GaussianCopula, n::Int64 = 1)
    L = cholesky(c.correlation).L

    Z = rand(Φ, n, size(c.correlation, 2))
    X = transpose(L * transpose(Z))

    u = cdf.(Φ, X)
end

function dimensions(c::GaussianCopula)
    size(c.correlation, 1)
end

function to_standard_normal_space(c::GaussianCopula, u)
    L = cholesky(c.correlation).L
    quantile.(Φ, u) * inv(L)'
end

function to_copula_space(c::GaussianCopula, s)
    L = cholesky(c.correlation).L
    cdf.(Φ, s * L')
end