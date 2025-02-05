"""
    AbstractPowerSpectralDensity

An abstract type representing a power spectral density model. Concrete types inheriting from this abstract type should implement specific models for power spectral density.
For the theoretical background of spectral densitz models see [Semi-empirical PSD functions](dynamics.md)

# Subtypes
- `CloughPenzien`
- `KanaiTajimi`

"""

abstract type AbstractPowerSpectralDensity end

"""
    CloughPenzien(w, S_0, ω_f, ζ_f, ω_g, ζ_g)

The `CloughPenzien` struct represents a power spectral density model with the following parameters:
- `w`: A vector of angular frequencies.
- `S_0`: A scaling factor.
- `ω_f`: Frequency parameter for the first oscillator.
- `ζ_f`: Damping ratio for the first oscillator.
- `ω_g`: Frequency parameter for the second oscillator.
- `ζ_g`: Damping ratio for the second oscillator.

The constructor calculates the power spectral density `p` based on the given parameters.

# Example
```julia
w = 0:0.1:10
S_0 = 1.0
ω_f = 2.0
ζ_f = 0.05
ω_g = 3.0
ζ_g = 0.1
cp = CloughPenzien(w, S_0, ω_f, ζ_f, ω_g, ζ_g)
```
"""
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

function evaluate(cp::CloughPenzien)
    return cp.p
end

function evaluate(kt::KanaiTajimi)
    return kt.S_0 .* (1 .+ 4 * kt.ζ^2 .* (kt.ω ./ kt.ω_0) .^ 2) ./
           ((1 .- (kt.ω ./ kt.ω_0) .^ 2) .^ 2 .+ 4 * kt.ζ^2 * (kt.ω ./ kt.ω_0) .^ 2)
end

struct ShinozukaDeodatis <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    σ::Real
    b::Real
    p::AbstractVector{<:Real}
end

function ShinozukaDeodatis(ω::AbstractVector{<:Real}, σ::Real, b::Real)
    p = 1 / 4 * σ^2 * b^3 .* ω .^ 2 .* exp.(-b * abs.(ω))
    return ShinozukaDeodatis(ω, σ, b, p)
end

function evaluate(sd::ShinozukaDeodatis)
    return sd.p
end

struct EmpiricalPSD <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    p::AbstractVector{<:Real}
end

function evaluate(ep::EmpiricalPSD)
    return ep.p
end
