# Power Spectral Density Functions

Construction and evaluation of different Power Spectral Density (PSD) functions. Correpsonding theory and literature can be found here [Semi-empirical-PSD-functions](@ref).

## Index

```@index
Pages = ["psd.md"]
```

## Types of PSD functions

```@docs
    CloughPenzien(ω::AbstractVector{<:Real}, S_0::Real, ω_f::Real, ζ_f::Real, ω_g::Real, ζ_g::Real)
    KanaiTajimi(ω::AbstractVector{<:Real}, S_0::Real, ω_0::Real, ζ::Real)
    ShinozukaDeodatis(ω::AbstractVector{<:Real}, σ::Real, b::Real)
    UncertaintyQuantification.EmpiricalPSD
```
