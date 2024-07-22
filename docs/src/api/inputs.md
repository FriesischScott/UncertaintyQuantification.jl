# Inputs

## Index

```@index
Pages = ["inputs.md"]
```

## Types

```@docs
Parameter
RandomVariable
Interval
ProbabilityBox
```

## Functions

```@docs
sample(rv::RandomVariable, n::Integer)
sample(inputs::Vector{<:UQInput}, n::Integer)

```
