# standard gaussian pdf
function ϕ(x::Real)
    return exp(-x^2 / 2) / sqrt(2 * π)
end
# 4th derivative of the standard gaussian pdf
function ϕ4(x::Real)
    return (exp(-x^2 / 2) * (x^4 - 6x^2 + 3)) / sqrt(2 * π)
end
# 6th derivative of the standard gaussian pdf
function ϕ6(x::Real)
    return (exp(-x^2 / 2) * (x^6 - 15x^4 + 45x^2 - 15)) / sqrt(2 * π)
end

function TD(b::Real, x::AbstractVector{<:Real})
    n = length(x)

    val = zero(eltype(x))

    @inbounds @simd for i in 1:(n - 1)
        for j in (i + 1):n
            val += @views 2 * ϕ6(b^(-1) * (x[i] - x[j]))
        end
    end

    val += n * ϕ6(zero(eltype(x)))

    return -1 * (n * (n - 1))^(-1) * b^(-7) * val
end

function SD(α::Real, x::AbstractVector)
    n = length(x)

    val = zero(eltype(x))

    @inbounds @simd for i in 1:(n - 1)
        for j in (i + 1):n
            val += @views 2 * ϕ4(α^(-1) * (x[i] - x[j]))
        end
    end

    val += n * ϕ4(zero(eltype(x)))

    return (n * (n - 1))^(-1) * α^(-5) * val
end

#=
S. J. Sheather, M. C. Jones, A Reliable Data-Based Bandwidth Selection Method for Kernel Density Estimation, Journal of the Royal Statistical Society: Series B (Methodological), Volume 53, Issue 3, July 1991, Pages 683–690, https://doi.org/10.1111/j.2517-6161.1991.tb01857.x
=#
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
