
# Polynomial Chaos Expansion {#Polynomial-Chaos-Expansion}

## Index {#Index}
- [`UncertaintyQuantification.WeightedApproximateFetekePoints`](#UncertaintyQuantification.WeightedApproximateFetekePoints)


## Types {#Types}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.WeightedApproximateFetekePoints' href='#UncertaintyQuantification.WeightedApproximateFetekePoints'><span class="jlbinding">UncertaintyQuantification.WeightedApproximateFetekePoints</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
WeightedApproximateFetekePoints(sim::AbstractMonteCarlo; fadd=10, fmult=2)
```


Struct for performing weighted approximate Feteke points (wafp) subsampling of a Monte-Carlo sampler for use in generating a  `PolynomialChaosExpansion`. Given a `PolynomialChaosBasis` of dimension `N`, and a Monte-Carlo sampler with `M` samples, generates a subsample of size `max(N,min(N+fadd,N+fmult,M))` biased towards maximizing the determinant of the Gramian typically requiring less  than `M` model evaluations. Follows procedure described in [[39](/references#burkEfficientSampling)]. 


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/pce/polynomialchaosexpansion.jl#L12-L19" target="_blank" rel="noreferrer">source</a></Badge>

</details>

