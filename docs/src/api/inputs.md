# Inputs

## Index

```@index
Pages = ["inputs.md"]
```

## Types

```@docs
Parameter
Interval
ProbabilityBox
GaussianMixtureModel
RandomVariable
IntervalVariable
EmpiricalDistribution
JointDistribution
```

## Functions

```@docs
sample(rv::RandomVariable, n::Integer=1)
sample(inputs::Vector{<:UQInput}, n::Integer=1)
sample(gmm::GaussianMixtureModel, n::Integer=1)
fit!(gmm::GaussianMixtureModel, data::DataFrame; max_iter::Integer=100, tol::Number=1e-4)
```
