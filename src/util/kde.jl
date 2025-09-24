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

function TD(b::Real, grid::AbstractVector{<:Real}, weights::AbstractVector{<:Real})
    n = sum(weights)
    nbins = length(grid)

    val = zero(eltype(grid))

    @inbounds @simd for i in 1:nbins-1
        for j in (i+1):nbins
                val += @views 2 * weights[i] * weights[j] * UncertaintyQuantification.ϕ6(b^(-1) * (grid[i] - grid[j]))
        end
    end

    val += sum(weights.^2) * UncertaintyQuantification.ϕ6(zero(eltype(grid)))

    return -1 * (n * (n - 1))^(-1) * b^(-7) * val
end

function SD(α::Real, grid::AbstractVector{<:Real}, weights::AbstractVector{<:Real})
    n = sum(weights)
    nbins = length(grid)

    val = zero(eltype(grid))

    @inbounds @simd for i in 1:(nbins - 1)
        for j in (i + 1):nbins
            val += @views 2 * weights[i] * weights[j] * UncertaintyQuantification.ϕ4(α^(-1) * (grid[i] - grid[j]))
        end
    end

    val += sum(weights.^2) * UncertaintyQuantification.ϕ4(zero(eltype(grid)))

    return (n * (n - 1))^(-1) * α^(-5) * val
end

#=
S. J. Sheather, M. C. Jones, A Reliable Data-Based Bandwidth Selection Method for Kernel Density Estimation, Journal of the Royal Statistical Society: Series B (Methodological), Volume 53, Issue 3, July 1991, Pages 683–690, https://doi.org/10.1111/j.2517-6161.1991.tb01857.x
=#
function sheather_jones_bandwidth(x::AbstractVector, nbins::Integer=0)
    λ = iqr(x)
    n = length(x)

    g, w = if nbins > 0
        linear_binning(x, nbins)
    else
        x, ones(n)
    end

    a = 0.920λ * n^(-1 / 7)
    b = 0.912λ * n^(-1 / 9)

    α₂ = 1.357 * (SD(a, g, w) / TD(b,g,w))^(1 / 7)

    function f(h)
        α₂_h = α₂ * h^(5 / 7)

        return ((1 / (2 * sqrt(π))) / (3 * SD(α₂_h, g, w)))^(1 / 5) * n^(-1 / 5) - h
    end

    return newtonraphson(1.0,f, 1e-3, 1e-5, 10^3)
end

function kde(h, x, X)
    return length(X)^(-1) * sum([h^(-1) * ϕ(h^(-1) * (x - xⱼ)) for xⱼ in X])
end
