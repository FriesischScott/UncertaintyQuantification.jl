```@meta
CurrentModule = UncertaintyQuantification
DocTestSetup = quote
    using UncertaintyQuantification
end
```

# Inputs

General functions operating on a collection of inputs defined as subtypes of `UQInput`.

## Index

```@index
Pages = ["inputs.md"]
Module = ["UncertaintyQuantification"]
```

## Functions
```@docs
sample(inputs::Array{<:UQInput}, n::Integer)
```