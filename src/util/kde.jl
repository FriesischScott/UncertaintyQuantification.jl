
# standard gaussian pdf
ϕ(x) = exp(-x^2 / 2) / sqrt(2 * π)
# 4th derivative of the standard gaussian pdf
ϕ4(x) = (exp(-x^2 / 2) * (x^4 - 6x^2 + 3)) / sqrt(2 * π)
# 6th derivative of the standard gaussian pdf
ϕ6(x) = (exp(-x^2 / 2) * (x^6 - 15x^4 + 45x^2 - 15)) / sqrt(2 * π)

function TD(b::Real, x::AbstractVector{<:Real})
    n = length(x)

    val = zero(eltype(x))

    for xᵢ in x, xⱼ in x
        val += ϕ6(b^(-1) * (xᵢ - xⱼ))
    end

    return -1 * (n * (n - 1))^(-1) * b^(-7) * val
end

function SD(α::Real, x::AbstractVector)
    n = length(x)

    val = zero(eltype(x))

    for xᵢ in x, xⱼ in x
        val += ϕ4(α^(-1) .* (xᵢ .- xⱼ))
    end

    return (n * (n - 1))^(-1) * α^(-5) * val
end

function sheather_jones_bandwidth(x::AbstractVector)
    λ = iqr(x)
    n = length(x)

    a = 0.920λ * n^(-1 / 7)
    b = 0.912λ * n^(-1 / 9)

    α₂ = 1.357 * (SD(a, x) / TD(b, x))^(1 / 7)

    function f(h, x)
        α₂_h = α₂ * h^(5 / 7)

        return ((1 / (2 * sqrt(π))) / (3 * SD(α₂_h, x)))^(1 / 5) * n^(-1 / 5) - h
    end

    return newtonraphson(1.0, h -> f(h, x), 1e-3, 1e-5, 10^3)
end

function kde(h, x, X)
    return length(X)^(-1) * sum([h^(-1) * ϕ(h^(-1) * (x - xⱼ)) for xⱼ in X])
end
