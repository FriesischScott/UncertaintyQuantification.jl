# Stochastic Processes (Spectral Representation)

Stochastic process generation based on the Spectral Representation Method which utilizes Power Spectral Density Functions. Correpsonding theory and literature can be found here [Stochastic-Process-Generation](@ref).

## Index

```@index
Pages = ["spectralrepresentation.md"]
```

## Types and Spectral Representation functions

```@docs
    SpectralRepresentation(psd::AbstractPowerSpectralDensity, time::AbstractVector{<:Real}, name::Symbol)
    sample(sr::SpectralRepresentation, n::Integer=1)
    evaluate(sr::SpectralRepresentation, ϕ::AbstractVector{<:Real})
    to_standard_normal_space!(sr::SpectralRepresentation, df::DataFrame)
    to_physical_space!(sr::SpectralRepresentation, df::DataFrame)
    dimensions(sr::SpectralRepresentation)
    names(sr::SpectralRepresentation)
```