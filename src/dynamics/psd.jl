
abstract type AbstractPowerSpectralDensity end

struct CloughPenzien <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    S_0::Real
    ω_f::Real
    ζ_f::Real
    ω_g::Real
    ζ_g::Real
    p::AbstractVector{<:Real}
end

function CloughPenzien(
    ω::AbstractVector{<:Real}, S_0::Real, ω_f::Real, ζ_f::Real, ω_g::Real, ζ_g::Real
)
    p =
        S_0 * ((ω .^ 4) ./ ((ω_f^2 .- ω .^ 2) .^ 2 .+ 4 * ζ_f^2 * ω_f^2 * ω .^ 2)) .* (
            (ω_g^4 .+ 4 * ζ_g^2 * ω_g^2 * ω .^ 2) ./
            ((ω_g^2 .- ω .^ 2) .^ 2 .+ (4 * ζ_g^2 * ω_g^2 * ω .^ 2))
        )

    return CloughPenzien(ω, S_0, ω_f, ζ_f, ω_g, ζ_g, p)
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

struct ShinozukaDeodatis <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    σ::Real
    b::Real
end

function evaluate(sd::ShinozukaDeodatis)
    return 1 / 4 * sd.σ^2 * sd.b^3 .* sd.ω .^ 2 .* exp.(-sd.b * abs.(sd.ω))
end

struct EmpiricalPSD <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    p::AbstractVector{<:Real}
end

function evaluate(ep::EmpiricalPSD)
    return ep.p
end
