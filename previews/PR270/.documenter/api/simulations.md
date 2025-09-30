
# Simulations {#Simulations}

Various Monte Carlo based simulations for a wide range of applications.

## Index {#Index}
- [`UncertaintyQuantification.RadialBasedImportanceSampling`](#UncertaintyQuantification.RadialBasedImportanceSampling)
- [`UncertaintyQuantification.SubSetInfinity`](#UncertaintyQuantification.SubSetInfinity)
- [`UncertaintyQuantification.SubSetInfinityAdaptive`](#UncertaintyQuantification.SubSetInfinityAdaptive)
- [`UncertaintyQuantification.SubSetSimulation`](#UncertaintyQuantification.SubSetSimulation)


## Types {#Types}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.RadialBasedImportanceSampling' href='#UncertaintyQuantification.RadialBasedImportanceSampling'><span class="jlbinding">UncertaintyQuantification.RadialBasedImportanceSampling</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
RadialBasedImportanceSampling(n::Integer, β::Real)
```


Used to perform radial-based importance sampling with `n` samples and reliability index `β`. If no `β` or `β=0.0` is passed, a [`FORM`](/api/reliability#UncertaintyQuantification.FORM) analysis will automatically be performed to estimate the reliability index.

**Examples**

`jldoctest julia> rbis = RadialBasedImportanceSampling(1000) RadialBasedImportanceSampling(10000, 0.0)``

**References**

[[15](/references#harbitzEfficientSamplingMethod1986)]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/simulations/radialbasedimportancesampling.jl#L1-L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.SubSetSimulation' href='#UncertaintyQuantification.SubSetSimulation'><span class="jlbinding">UncertaintyQuantification.SubSetSimulation</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SubSetSimulation(n::Integer, target::Float64, levels::Integer, proposal::UnivariateDistribution)
```


Defines the properties of a Subset simulation where `n` is the number of initial samples, `target` is the target probability of failure at each level, `levels` is the maximum number of levels and `proposal` is the proposal distribution for the markov chain monte carlo.

**Examples**

```julia
julia> SubSetSimulation(100, 0.1, 10, Uniform(-0.2, 0.2))
SubSetSimulation(100, 0.1, 10, Uniform{Float64}(a=-0.2, b=0.2))
```


**References**

[[18](/references#auEstimationSmallFailure2001)]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/simulations/subset.jl#L3-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.SubSetInfinity' href='#UncertaintyQuantification.SubSetInfinity'><span class="jlbinding">UncertaintyQuantification.SubSetInfinity</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SubSetInfinity(n::Integer, target::Float64, levels::Integer, s::Real)
```


Defines the properties of a Subset-∞ simulation where `n` is the number of initial samples, `target` is the target probability of failure at each level, `levels` is the maximum number of levels and `s` is the standard deviation for the proposal samples.

**Examples**

```julia
julia> SubSetInfinity(100, 0.1, 10, 0.5)
SubSetInfinity(100, 0.1, 10, 0.5)
```


**References**

[[19](/references#auRareEventSimulation2016)]

[[36](/references#patelliEfficientMonteCarlo2015)]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/simulations/subset.jl#L43-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.SubSetInfinityAdaptive' href='#UncertaintyQuantification.SubSetInfinityAdaptive'><span class="jlbinding">UncertaintyQuantification.SubSetInfinityAdaptive</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real, s::Real)
```


Implementation of: Papaioannou, Iason, et al. &quot;MCMC algorithms for subset simulation.&quot; Probabilistic Engineering Mechanics 41 (2015): 89-103

Defines the properties of a Subset-∞ adaptive where `n` are the number of samples per level, `target` is the target probability of failure at each level, `levels` is the maximum number of levels, `λ` (λ = 1 recommended) is the initial scaling parameter, and `Na` is the number simulations that will be run before `λ` is updated. Note that Na must be a multiple of n * target: `mod(ceil(n * target), Na) == 0)`. The initial variance of the proposal distribution is `s`.

Idea behind this algorithm is to adaptively select the correlation parameter of `s` at each intermediate level, by simulating a subset Na of the chains (which must be choosen without replacement at random) and modifying the acceptance rate towards the optimal αstar = 0.44

**Constructors**
- `SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer)`   (default: λ = s = 1)
  
- `SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real)` (λ = s)
  
- `SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real, s::Real)`
  

**Note**

The following constructors will run the same number of samples, but SubSetInfinityAdaptive will update `s` after each chain:
- `SubSetInfinityAdaptive(400, 0.1, 10, 40)`
  
- `SubSetInfinity(400, 0.1, 10, 0.5)`
  

**Examples**

```julia
julia> SubSetInfinityAdaptive(200, 0.1, 10, 2)
SubSetInfinityAdaptive(200, 0.1, 10, 2, 1, 1)
```


**References**

[[37](/references#papaioannou2015mcmc)]

[[38](/references#chan2022adaptive)]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/simulations/subset.jl#L92-L132" target="_blank" rel="noreferrer">source</a></Badge>

</details>

