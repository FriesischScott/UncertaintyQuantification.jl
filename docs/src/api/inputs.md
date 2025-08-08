# Inputs

## Index

```@index
Pages = ["inputs.md"]
```

## Types

```@docs
Parameter
RandomVariable
EmpiricalDistribution
Interval
ProbabilityBox
GaussianMixtureModel
```

## Functions

```@docs
sample(rv::RandomVariable, n::Integer=1)
sample(inputs::Vector{<:UQInput}, n::Integer=1)
sample(gmm::GaussianMixtureModel, n::Integer=1)
fit_gaussian_mixture(number_components::Integer, data::Matrix; maximum_iterations::Integer=100, tolerance::Number=1e-4)
```
