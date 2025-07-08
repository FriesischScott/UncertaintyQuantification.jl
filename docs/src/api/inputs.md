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
IntervalVariable
ProbabilityBox
```

## Functions

```@docs
sample(rv::RandomVariable, n::Integer=1)
sample(inputs::Vector{<:UQInput}, n::Integer=1)

```
