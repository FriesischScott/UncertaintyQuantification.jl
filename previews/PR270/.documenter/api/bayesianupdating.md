
# Bayesian Updating {#Bayesian-Updating}

Methods for Bayesian updating.

## Index {#Index}
- [`UncertaintyQuantification.AbstractBayesianMethod`](#UncertaintyQuantification.AbstractBayesianMethod)
- [`UncertaintyQuantification.AbstractBayesianPointEstimate`](#UncertaintyQuantification.AbstractBayesianPointEstimate)
- [`UncertaintyQuantification.MaximumAPosterioriBayesian`](#UncertaintyQuantification.MaximumAPosterioriBayesian)
- [`UncertaintyQuantification.MaximumLikelihoodBayesian`](#UncertaintyQuantification.MaximumLikelihoodBayesian)
- [`UncertaintyQuantification.SingleComponentMetropolisHastings`](#UncertaintyQuantification.SingleComponentMetropolisHastings)
- [`UncertaintyQuantification.TransitionalMarkovChainMonteCarlo`](#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo)
- [`UncertaintyQuantification.bayesianupdating`](#UncertaintyQuantification.bayesianupdating)


## Types {#Types}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.AbstractBayesianMethod' href='#UncertaintyQuantification.AbstractBayesianMethod'><span class="jlbinding">UncertaintyQuantification.AbstractBayesianMethod</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
AbstractBayesianMethod
```


Subtypes are used to dispatch to the different MCMC methods in [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating).

Subtypes are:
- [`SingleComponentMetropolisHastings`](/api/bayesianupdating#UncertaintyQuantification.SingleComponentMetropolisHastings)
  
- [`TransitionalMarkovChainMonteCarlo`](/api/bayesianupdating#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/UncertaintyQuantification.jl#L46-L55" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.SingleComponentMetropolisHastings' href='#UncertaintyQuantification.SingleComponentMetropolisHastings'><span class="jlbinding">UncertaintyQuantification.SingleComponentMetropolisHastings</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SingleComponentMetropolisHastings(proposal, x0, n, burnin, islog)
```


Passed to [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating) to run the single-component Metropolis-Hastings algorithm starting from `x0` with  univariate proposal distibution `proposal`. Will generate `n` samples _after_ performing `burnin` steps of the Markov chain and discarding the samples. The flag `islog` specifies whether the prior and likelihood functions passed to the  [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating) method are already  given as logarithms.

Alternative constructor

```julia
    SingleComponentMetropolisHastings(proposal, x0, n, burnin)  # `islog` = true
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/modelupdating/bayesianupdating.jl#L1-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.TransitionalMarkovChainMonteCarlo' href='#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo'><span class="jlbinding">UncertaintyQuantification.TransitionalMarkovChainMonteCarlo</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TransitionalMarkovChainMonteCarlo(prior, n, burnin, β, islog)

Passed to [`bayesianupdating`](@ref) to run thetransitional Markov chain Monte Carlo algorithm  with [`RandomVariable'](@ref) vector `prior`. At each transitional level, one sample will be generated from `n` independent Markov chains after `burnin` steps have been discarded. The flag `islog` specifies whether the prior and likelihood functions passed to the  [`bayesianupdating`](@ref) method are already  given as logarithms.
```


Alternative constructors

```julia
    TransitionalMarkovChainMonteCarlo(prior, n, burnin, β)  # `islog` = true
     TransitionalMarkovChainMonteCarlo(prior, n, burnin)    # `β` = 0.2,  `islog` = true
```


**References**

[[30](/references#chingTransitionalMarkovChain2007)]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/modelupdating/bayesianupdating.jl#L134-L150" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.AbstractBayesianPointEstimate' href='#UncertaintyQuantification.AbstractBayesianPointEstimate'><span class="jlbinding">UncertaintyQuantification.AbstractBayesianPointEstimate</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
AbstractBayesianPointEstimate
```


Subtypes are used to dispatch to the different point estimation methods in [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating).

Subtypes are:
- [`MaximumAPosterioriBayesian`](/api/bayesianupdating#UncertaintyQuantification.MaximumAPosterioriBayesian)
  
- [`MaximumLikelihoodBayesian`](/api/bayesianupdating#UncertaintyQuantification.MaximumLikelihoodBayesian)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/UncertaintyQuantification.jl#L58-L67" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.MaximumAPosterioriBayesian' href='#UncertaintyQuantification.MaximumAPosterioriBayesian'><span class="jlbinding">UncertaintyQuantification.MaximumAPosterioriBayesian</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MaximumAPosterioriBayesian(prior, optimmethod, x0; islog, lowerbounds, upperbounds)
```


Passed to [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating) to estimate one or more maxima of the posterior distribution starting from `x0`. The optimization uses the method specified in `optimmethod`. Will calculate one estimation per point in x0. The flag `islog` specifies whether the prior and likelihood functions passed to the  [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating) method are already  given as logarithms. `lowerbounds` and `upperbounds` specify optimization intervals.

Alternative constructors

```julia
    MaximumAPosterioriBayesian(prior, optimmethod, x0; islog) # `lowerbounds` = [-Inf], # `upperbounds` = [Inf]
    MaximumAPosterioriBayesian(prior, optimmethod, x0)  # `islog` = true
```


See also [`MaximumLikelihoodBayesian`](/api/bayesianupdating#UncertaintyQuantification.MaximumLikelihoodBayesian), [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating),  [`TransitionalMarkovChainMonteCarlo`](/api/bayesianupdating#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/modelupdating/bayesianMAP.jl#L1-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.MaximumLikelihoodBayesian' href='#UncertaintyQuantification.MaximumLikelihoodBayesian'><span class="jlbinding">UncertaintyQuantification.MaximumLikelihoodBayesian</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MaximumLikelihoodBayesian(prior, optimmethod, x0; islog, lowerbounds, upperbounds)
```


Passed to [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating) to estimate one or more maxima of the likelihood starting from `x0`. The optimization uses the method specified in `optimmethod`. Will calculate one estimation per point in x0. The flag `islog` specifies whether the prior and likelihood functions passed to the  [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating) method are already  given as logarithms. `lowerbounds` and `upperbounds` specify optimization intervals.

Alternative constructors

```julia
    MaximumLikelihoodBayesian(prior, optimmethod, x0; islog) # `lowerbounds` = [-Inf], # `upperbounds` = [Inf]
    MaximumLikelihoodBayesian(prior, optimmethod, x0)  # `islog` = true
```


**Notes**

The method uses `prior` only as information on which parameters are supposed to be optimized. The prior itself does not influence the result of the maximum likelihood estimate and can be given as a dummy distribution. For example, if two parameters `a` and `b` are supposed to be optimized, the prior could look like this

```julia
    prior = RandomVariable.(Uniform(0,1), [:a, :b])
```


See also [`MaximumAPosterioriBayesian`](/api/bayesianupdating#UncertaintyQuantification.MaximumAPosterioriBayesian), [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating),  [`TransitionalMarkovChainMonteCarlo`](/api/bayesianupdating#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/modelupdating/bayesianMAP.jl#L52-L70" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Methods {#Methods}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.bayesianupdating' href='#UncertaintyQuantification.bayesianupdating'><span class="jlbinding">UncertaintyQuantification.bayesianupdating</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
bayesianupdating(likelihood, models, pointestimate; prior)
```


Perform bayesian updating using the given `likelihood`, `models`  and any point estimation method [`AbstractBayesianPointEstimate`](/api/bayesianupdating#UncertaintyQuantification.AbstractBayesianPointEstimate).

**Notes**

Method can be called with an empty Vector of models, i.e.

```
bayesianupdating(likelihood, [], pointestimate)
```


If `prior` is not given, the method will construct a prior distribution from the prior specified in `AbstractBayesianPointEstimate.prior`.

`likelihood` is a Julia function which must be defined in terms of a `DataFrame` of samples, and must evaluate the likelihood for each row of the `DataFrame`

For example, a loglikelihood based on normal distribution using &#39;Data&#39;:

```julia
likelihood(df) = [sum(logpdf.(Normal.(df_i.x, 1), Data)) for df_i in eachrow(df)]
```


If a model evaluation is required to evaluate the likelihood, a vector of `UQModel`s must be passed to `bayesianupdating`. For example if the variable `x` above is the output of a numerical model.

For a general overview of the function, see [`bayesianupdating`](/api/bayesianupdating#UncertaintyQuantification.bayesianupdating).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/modelupdating/bayesianMAP.jl#L112-L136" target="_blank" rel="noreferrer">source</a></Badge>



```julia
bayesianupdating(prior, likelihood, models, mcmc)
```


Perform bayesian updating using the given `prior`, `likelihood`, `models`  and any MCMC sampler [`AbstractBayesianMethod`](/api/bayesianupdating#UncertaintyQuantification.AbstractBayesianMethod).

Alternatively the method can be called without `models`.

```
bayesianupdating(prior, likelihood, mcmc)
```


When using [`TransitionalMarkovChainMonteCarlo`](/api/bayesianupdating#UncertaintyQuantification.TransitionalMarkovChainMonteCarlo) the `prior` can automatically be constructed.

```
bayesinupdating(likelihood, models, tmcmc)
bayesianupdating(likelihood, tmcmc)
```


**Notes**

`likelihood` is a Julia function which must be defined in terms of a `DataFrame` of samples, and must evaluate the likelihood for each row of the `DataFrame`

For example, a loglikelihood based on normal distribution using &#39;Data&#39;:

```julia
likelihood(df) = [sum(logpdf.(Normal.(df_i.x, 1), Data)) for df_i in eachrow(df)]
```


If a model evaluation is required to evaluate the likelihood, a vector of `UQModel`s must be passed to `bayesianupdating`. For example if the variable `x` above is the output of a numerical model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/modelupdating/bayesianupdating.jl#L34-L61" target="_blank" rel="noreferrer">source</a></Badge>

</details>

