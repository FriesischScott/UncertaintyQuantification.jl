# Reliability

## Index

```@index
Pages = ["reliability.md"]
```

## Types

```@docs
FORM
```

## Methods

```@docs
probability_of_failure(models::Union{Vector{<:UQModel},UQModel},performance::Function,inputs::Union{Vector{<:UQInput},UQInput},sim::FORM)
probability_of_failure(models::Union{Vector{<:UQModel},UQModel},performance::Function,inputs::Union{Vector{<:UQInput},UQInput},sim::AbstractMonteCarlo)
```
