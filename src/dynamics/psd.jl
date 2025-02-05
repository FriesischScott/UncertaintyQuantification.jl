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

"""
    CloughPenzien(ω::AbstractVector{<:Real}, S_0::Real, ω_f::Real, ζ_f::Real, ω_g::Real, ζ_g::Real)

Constructs a `CloughPenzien` instance representing a power spectral density model with the given parameters.

# Arguments
- `ω::AbstractVector{<:Real}`: A vector of angular frequencies.
- `S_0::Real`: A scaling factor.
- `ω_f::Real`: Frequency parameter for the first oscillator.
- `ζ_f::Real`: Damping ratio for the first oscillator.
- `ω_g::Real`: Frequency parameter for the second oscillator.
- `ζ_g::Real`: Damping ratio for the second oscillator.

# Returns
A `CloughPenzien` instance with the calculated power spectral density `p`.

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

"""
    KanaiTajimi(ω::AbstractVector{<:Real}, S_0::Real, ω_0::Real, ζ::Real)

The `KanaiTajimi` struct represents a power spectral density model with the following parameters:
- `ω`: A vector of angular frequencies.
- `S_0`: A scaling factor.
- `ω_0`: Natural frequency of the oscillator.
- `ζ`: Damping ratio of the oscillator.

# Example
```julia
w = 0:0.1:10
S_0 = 1.0
ω_0 = 2.0
ζ = 0.05
kt = KanaiTajimi(w, S_0, ω_0, ζ)
```
"""
struct KanaiTajimi <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    S_0::Real
    ω_0::Real
    ζ::Real
end

"""
    evaluate(cp::CloughPenzien) -> AbstractVector{<:Real}

Evaluates the power spectral density for a given `CloughPenzien` instance.

# Arguments
- `cp::CloughPenzien`: An instance of the `CloughPenzien` struct.

# Returns
A vector of real numbers representing the values of the power spectral density function.

# Example
```julia
w = 0:0.1:10
S_0 = 1.0
ω_f = 2.0
ζ_f = 0.05
ω_g = 3.0
ζ_g = 0.1
cp = CloughPenzien(w, S_0, ω_f, ζ_f, ω_g, ζ_g)
psd_values = evaluate(cp)
```
"""
function evaluate(cp::CloughPenzien)
    return cp.p
end

"""
    evaluate(kt::KanaiTajimi) -> AbstractVector{<:Real}

Evaluates the power spectral density for a given `KanaiTajimi` instance.

# Arguments
- `kt::KanaiTajimi`: An instance of the `KanaiTajimi` struct.

# Returns
A vector of real numbers representing the power spectral density.

# Example
```julia
w = 0:0.1:10
S_0 = 1.0
ω_0 = 2.0
ζ = 0.05
kt = KanaiTajimi(w, S_0, ω_0, ζ)
psd_values = evaluate(kt)
```
"""
function evaluate(kt::KanaiTajimi)
    return kt.S_0 .* (1 .+ 4 * kt.ζ^2 .* (kt.ω ./ kt.ω_0) .^ 2) ./
           ((1 .- (kt.ω ./ kt.ω_0) .^ 2) .^ 2 .+ 4 * kt.ζ^2 * (kt.ω ./ kt.ω_0) .^ 2)
end

"""
    ShinozukaDeodatis(ω::AbstractVector{<:Real}, σ::Real, b::Real)

The `ShinozukaDeodatis` struct represents a power spectral density model with the following parameters:
- `ω`: A vector of angular frequencies.
- `σ`: A scaling factor.
- `b`: A parameter related to the model.
- `p`: A vector representing the power spectral density, which will be computed based on the given parameters.

# Example
```julia
w = 0:0.1:10
σ = 1.0
b = 0.5
sd = ShinozukaDeodatis(w, σ, b)
```
"""
struct ShinozukaDeodatis <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    σ::Real
    b::Real
    p::AbstractVector{<:Real}
end

"""
    ShinozukaDeodatis(ω::AbstractVector{<:Real}, σ::Real, b::Real)

Constructs a `ShinozukaDeodatis` instance representing a power spectral density model with the given parameters.

# Arguments
- `ω::AbstractVector{<:Real}`: A vector of angular frequencies.
- `σ::Real`: A scaling factor.
- `b::Real`: A parameter related to the model.

# Returns
A `ShinozukaDeodatis` instance with the calculated power spectral density `p`.

# Example
```julia
w = 0:0.1:10
σ = 1.0
b = 0.5
sd = ShinozukaDeodatis(w, σ, b)
```
"""
function ShinozukaDeodatis(ω::AbstractVector{<:Real}, σ::Real, b::Real)
    p = 1 / 4 * σ^2 * b^3 .* ω .^ 2 .* exp.(-b * abs.(ω))
    return ShinozukaDeodatis(ω, σ, b, p)
end

"""
    evaluate(sd::ShinozukaDeodatis) -> AbstractVector{<:Real}

Evaluates the power spectral density for a given `ShinozukaDeodatis` instance.

# Arguments
- `sd::ShinozukaDeodatis`: An instance of the `ShinozukaDeodatis` struct.

# Returns
A vector of real numbers representing the power spectral density.

# Example
```julia
w = 0:0.1:10
σ = 1.0
b = 0.5
sd = ShinozukaDeodatis(w, σ, b)
psd_values = evaluate(sd)
```
"""
function evaluate(sd::ShinozukaDeodatis)
    return sd.p
end

"""
    EmpiricalPSD(ω::AbstractVector{<:Real}, p::AbstractVector{<:Real})

The `EmpiricalPSD` struct represents an empirical power spectral density model with the following parameters:
- `ω`: A vector of angular frequencies.
- `p`: A vector of power spectral density values corresponding to the frequencies in `ω`.

# Example
```julia
w = 0:0.1:10
p_values = rand(length(w))  # Example empirical PSD values
emp_psd = EmpiricalPSD(w, p_values)
```
"""
struct EmpiricalPSD <: AbstractPowerSpectralDensity
    ω::AbstractVector{<:Real}
    p::AbstractVector{<:Real}
end

"""
    evaluate(ep::EmpiricalPSD) -> AbstractVector{<:Real}

Evaluates the power spectral density for a given `EmpiricalPSD` instance.

# Arguments
- `ep::EmpiricalPSD`: An instance of the `EmpiricalPSD` struct.

# Returns
A vector of real numbers representing the power spectral density.

# Example
```julia
w = 0:0.1:10
p_values = rand(length(w))  # Example empirical PSD values
emp_psd = EmpiricalPSD(w, p_values)
psd_values = evaluate(emp_psd)
```
"""
function evaluate(ep::EmpiricalPSD)
    return ep.p
end