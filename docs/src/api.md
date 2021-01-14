# API

```@meta
CurrentModule = UncertaintyQuantification
DocTestSetup = quote
    using UncertaintyQuantification
end
```

## Inputs

```@docs
Parameter
RandomVariable
sample(inputs::Array{<:UQInput}, n::Integer)
sample(rv::RandomVariable, n::Integer)

```

## Index

```@index
Pages = ["api.md"]
Module = ["UncertaintyQuantification"]
```