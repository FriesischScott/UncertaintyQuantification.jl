
# Reliability {#Reliability}

## Index {#Index}
- [`UncertaintyQuantification.DoubleLoop`](#UncertaintyQuantification.DoubleLoop)
- [`UncertaintyQuantification.FORM`](#UncertaintyQuantification.FORM)
- [`UncertaintyQuantification.RandomSlicing`](#UncertaintyQuantification.RandomSlicing)
- [`UncertaintyQuantification.probability_of_failure`](#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel,%20Vector{<:UQModel}},%20Function,%20Union{UQInput,%20Vector{<:UQInput}},%20AbstractMonteCarlo})
- [`UncertaintyQuantification.probability_of_failure`](#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel,%20Vector{<:UQModel}},%20Function,%20Union{UQInput,%20Vector{<:UQInput}},%20FORM})
- [`UncertaintyQuantification.probability_of_failure`](#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel,%20Vector{<:UQModel}},%20Function,%20Union{UQInput,%20Vector{<:UQInput}},%20RandomSlicing})
- [`UncertaintyQuantification.probability_of_failure`](#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel,%20Vector{<:UQModel}},%20Function,%20Union{UQInput,%20Vector{<:UQInput}},%20DoubleLoop})


## Types {#Types}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.FORM' href='#UncertaintyQuantification.FORM'><span class="jlbinding">UncertaintyQuantification.FORM</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
FORM(n::Integer=10,tol::Real=1e-3,fdm::FiniteDifferencesMethod=CentralFiniteDifferences(3))
```


used to perform the first order reliability method using the HLRF algorithm with `n` iterations and tolerance `tol`. Gradients are estimated through `fdm`.

**References**

[[13](/references#rackwitzStructuralReliability1978)]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/reliability/form.jl#L1-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.DoubleLoop' href='#UncertaintyQuantification.DoubleLoop'><span class="jlbinding">UncertaintyQuantification.DoubleLoop</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DoubleLoop(lb::AbstractSimulation, ub::AbstractSimulation)
```


Used to estimate imprecise reliability with the _double loop_ Monte Carlo method. 

Wraps two simulation objects — one for lower-bound (`lb`) and one for upper-bound (`ub`).

The two simulations can differ in simulation type, complexity, or accuracy settings, since estimating the lower bound often requires more simulation effort.

This approach runs an optimisation loop over interval parameters (outer loop) and computes reliability bounds in an inner loop using the `lb` and `ub` simulation methods.

Use `DoubleLoop(sim::AbstractSimulation)` for creating a `DoubleLoop` with same simulation method for both bounds.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/reliability/probabilityoffailure_imprecise.jl#L2-L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.RandomSlicing' href='#UncertaintyQuantification.RandomSlicing'><span class="jlbinding">UncertaintyQuantification.RandomSlicing</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
RandomSlicing(lb::AbstractSimulation, ub::AbstractSimulation)
```


Used to estimate imprecise reliability with _random slicing_ Monte Carlo method, sometimes known as interval Monte Carlo.

Wraps two simulation objects — one for lower-bound (`lb`) and one for upper-bound (`ub`). 

The two simulations can differ in simulation type, complexity, or accuracy settings, since estimating the lower bound often requires more simulation effort.

In this approach, the `lb` and `ub` simulation methods generate random intervals from the imprecise variables. These intervals are then propagated through the model via optimisation-based interval propagation, yielding lower and upper bounds on the reliability estimate.

Use `RandomSlicing(sim::AbstractSimulation)` for creating a `RandomSlicing` with same simulation method for both bounds.

**References**

[[21](/references#alvarez2018estimation)]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/reliability/probabilityoffailure_imprecise.jl#L30-L46" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Methods {#Methods}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, FORM}' href='#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, FORM}'><span class="jlbinding">UncertaintyQuantification.probability_of_failure</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
probability_of_failure(models::Union{Vector{<:UQModel},UQModel},performance::Function),inputs::Union{Vector{<:UQInput},UQInput},sim::FORM)
```


Perform a reliability analysis using the first order reliability method (FORM), see [`FORM`](/api/reliability#UncertaintyQuantification.FORM). Returns the estimated probability of failure `pf`, the reliability index `β` and the design point `dp`.

**Examples**

```
pf, β, dp = probability_of_failure(model, performance, inputs, sim)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/reliability/form.jl#L24-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, AbstractMonteCarlo}' href='#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, AbstractMonteCarlo}'><span class="jlbinding">UncertaintyQuantification.probability_of_failure</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
probability_of_failure(models::Union{Vector{<:UQModel},UQModel},performance::Function),inputs::Union{Vector{<:UQInput},UQInput},sim::AbstractMonteCarlo)
```


Perform a reliability analysis with a standard Monte Carlo simulation. Returns the estimated probability of failure `pf`, the standard deviation `σ` and the `DataFrame` containing the evaluated `samples`. The simulation `sim` can be any instance of `AbstractMonteCarlo`.

**Examples**

```
pf, σ, samples = probability_of_failure(model, performance, inputs, sim)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/reliability/probabilityoffailure.jl#L1-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, DoubleLoop}' href='#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, DoubleLoop}'><span class="jlbinding">UncertaintyQuantification.probability_of_failure</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
probability_of_failure(
    models::Union{Vector{<:UQModel}, UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput}, UQInput},
    dl::DoubleLoop
)
```


Perform an **imprecise reliability analysis** using the _double loop_ Monte Carlo method.

The inputs must include at least one imprecise variable.

**Returns**
- **`pf_bounds`**: An [`Interval`](/api/inputs#UncertaintyQuantification.Interval) giving the lower and upper bounds on the probability of failure.  
  
- **`result_lb`**: The outputs of the reliability simulation that achieved the lower bound.  
  
- **`result_ub`**: The outputs of the reliability simulation that achieved the upper bound. 
  

If the lower and upper bounds are equal, returns only the scalar probability of failure.

See [`DoubleLoop`](/api/reliability#UncertaintyQuantification.DoubleLoop) for details of the random slicing configuration.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/reliability/probabilityoffailure_imprecise.jl#L62-L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, RandomSlicing}' href='#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel, Vector{<:UQModel}}, Function, Union{UQInput, Vector{<:UQInput}}, RandomSlicing}'><span class="jlbinding">UncertaintyQuantification.probability_of_failure</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
probability_of_failure(
    models::Union{Vector{<:UQModel}, UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput}, UQInput},
    rs::RandomSlicing
)
```


Perform an **imprecise reliability analysis** using the _random slicing_ Monte Carlo method

The inputs must include at least one imprecise variable.  

**Returns**
- **`pf_bounds`**: An [`Interval`](/api/inputs#UncertaintyQuantification.Interval) giving the lower and upper bounds on the probability of failure.  
  
- **`result_lb`**: The outputs of the reliability simulation that achieved the lower bound.  
  
- **`result_ub`**: The outputs of the reliability simulation that achieved the upper bound.
  

See [`RandomSlicing`](/api/reliability#UncertaintyQuantification.RandomSlicing) for details of the random slicing configuration.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/reliability/probabilityoffailure_imprecise.jl#L175-L193" target="_blank" rel="noreferrer">source</a></Badge>

</details>

