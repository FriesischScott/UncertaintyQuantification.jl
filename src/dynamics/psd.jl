
abstract type AbstractPowerSpectralDensity end

struct CloughPenzien <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    S_0::Real
    ω_f::Real
    ζ_f::Real
    ω_g::Real
    ζ_g::Real
end

function evaluate(cp::CloughPenzien)
    return cp.S_0 * (
        (cp.ω .^ 4) ./ ((cp.ω_f^2 .- cp.ω .^ 2) .^ 2 .+ 4 * cp.ζ_f^2 * cp.ω_f^2 * cp.ω .^ 2)
    ) .* (
        (cp.ω_g^4 .+ 4 * cp.ζ_g^2 * cp.ω_g^2 * cp.ω .^ 2) ./
        ((cp.ω_g^2 .- cp.ω .^ 2) .^ 2 .+ (4 * cp.ζ_g^2 * cp.ω_g^2 * cp.ω .^ 2))
    )
end

struct KanaiTajimi <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    S_0::Real
    ω_0::Real
    ζ::Real
end

function evaluate(kt::KanaiTajimi)
    return kt.S_0 .* (1 .+ 4 * kt.ζ^2 .* (kt.ω ./ kt.ω_0) .^ 2) ./
           ((1 .- (kt.ω ./ kt.ω_0) .^ 2) .^ 2 .+ 4 * kt.ζ^2 * (kt.ω ./ kt.ω_0) .^ 2)
end
