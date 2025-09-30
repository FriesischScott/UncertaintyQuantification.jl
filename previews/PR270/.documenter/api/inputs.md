
# Inputs {#Inputs}

## Index {#Index}
- [`UncertaintyQuantification.EmpiricalDistribution`](#UncertaintyQuantification.EmpiricalDistribution)
- [`UncertaintyQuantification.Interval`](#UncertaintyQuantification.Interval)
- [`UncertaintyQuantification.IntervalVariable`](#UncertaintyQuantification.IntervalVariable)
- [`UncertaintyQuantification.JointDistribution`](#UncertaintyQuantification.JointDistribution)
- [`UncertaintyQuantification.Parameter`](#UncertaintyQuantification.Parameter)
- [`UncertaintyQuantification.ProbabilityBox`](#UncertaintyQuantification.ProbabilityBox)
- [`UncertaintyQuantification.RandomVariable`](#UncertaintyQuantification.RandomVariable)
- [`UncertaintyQuantification.GaussianMixtureModel`](#UncertaintyQuantification.GaussianMixtureModel)
- [`UncertaintyQuantification.sample`](#UncertaintyQuantification.sample)
- [`UncertaintyQuantification.sample`](#UncertaintyQuantification.sample)


## Types {#Types}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.Parameter' href='#UncertaintyQuantification.Parameter'><span class="jlbinding">UncertaintyQuantification.Parameter</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Parameter(value::Real, name::Symbol)
```


Defines a parameter value (scalar), with an input value and a name.

**Examples**

```julia
julia> Parameter(3.14, :π)
Parameter(3.14, :π)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/parameter.jl#L1-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.Interval' href='#UncertaintyQuantification.Interval'><span class="jlbinding">UncertaintyQuantification.Interval</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Interval(lb::Real, ub::Real)
```


Represents a closed numeric interval with a lower bound `lb` and an upper bound `ub`.

`Interval` is a data type primarily used for constructing probability boxes (p-boxes) and other uncertainty representations. It is **not** intended for direct use in simulations for that, see [`IntervalVariable`](/api/inputs#UncertaintyQuantification.IntervalVariable).

**Fields**
- `lb::Real`: Lower bound of the interval.
  
- `ub::Real`: Upper bound of the interval.
  

**Examples**

```julia
julia> Interval(0.10, 0.14)
[0.1, 0.14]
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/imprecise/interval.jl#L1-L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.ProbabilityBox' href='#UncertaintyQuantification.ProbabilityBox'><span class="jlbinding">UncertaintyQuantification.ProbabilityBox</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ProbabilityBox{T}(parameters::Dict{Symbol, Union{Real, Interval}}, lb::Real, ub::Real)
```


Represents a (optionally truncated) probability box (p-box) for a univariate distribution `T`, where parameters may be specified as precise values (`Real`) or intervals ([`Interval`](/api/inputs#UncertaintyQuantification.Interval)). The support of the distribution is bounded by `lb` (lower bound) and `ub` (upper bound).

To use the `ProbabilityBox` in an analysis it has to be wrapped in a [`RandomVariable`](/api/inputs#UncertaintyQuantification.RandomVariable).

**Fields**
- `parameters::Dict{Symbol, Union{Real, Interval}}`: Dictionary mapping parameter names (as symbols) to their values or intervals.
  
- `lb::Real`: Lower bound of the distribution&#39;s support.
  
- `ub::Real`: Upper bound of the distribution&#39;s support.
  

**Constructors**
- `ProbabilityBox{T}(parameters::Dict{Symbol, Union{Real, Interval}}, lb::Real, ub::Real)`: Specify all parameters and support bounds explicitly.
  
- `ProbabilityBox{T}(parameters::Dict{Symbol, Union{Real, Interval}})`: Support bounds are inferred from the distribution type `T`.
  
- `ProbabilityBox{T}(parameter::Interval)`: For univariate distributions with a single parameter, construct from a single interval.
  

For convenience the first two constructors can also be called with a `Vector{Pair{Symbol,Union{Real,Interval}}}` to automatically create the `Dict`.

**Examples**

```julia
julia> ProbabilityBox{Normal}(Dict(:μ => Interval(0, 1), :σ => Interval(0.1, 1)), 0.0, Inf)
ProbabilityBox{Normal}(Dict{Symbol, Union{Real, Interval}}(:μ => [0, 1], :σ => [0.1, 1]), 0.0, Inf)

julia> ProbabilityBox{Normal}(Dict(:μ => Interval(0, 1), :σ => Interval(0.1, 1)))
ProbabilityBox{Normal}(Dict{Symbol, Union{Real, Interval}}(:μ => [0, 1], :σ => [0.1, 1]), -Inf, Inf)

julia> ProbabilityBox{Exponential}(Interval(0.1, 0.5))
ProbabilityBox{Exponential}(Dict{Symbol, Union{Real, Interval}}(:θ => [0.1, 0.5]), 0.0, Inf)

julia> ProbabilityBox{Normal}([:μ => Interval(0, 1), :σ => Interval(0.1, 1)])
ProbabilityBox{Normal}(Dict{Symbol, Union{Real, Interval}}(:μ => [0, 1], :σ => [0.1, 1]), -Inf, Inf)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/imprecise/p-box.jl#L1-L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.RandomVariable' href='#UncertaintyQuantification.RandomVariable'><span class="jlbinding">UncertaintyQuantification.RandomVariable</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
RandomVariable(dist::UnivariateDistribution, name::Symbol)
```


Defines a random variable, with a univariate distribution from Distributions.jl and a name.

**Examples**

```julia
julia> RandomVariable(Normal(), :x)
RandomVariable{Normal{Float64}}(Normal{Float64}(μ=0.0, σ=1.0), :x)

julia> RandomVariable(Exponential(1), :x)
RandomVariable{Exponential{Float64}}(Exponential{Float64}(θ=1.0), :x)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/randomvariables/randomvariable.jl#L1-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.IntervalVariable' href='#UncertaintyQuantification.IntervalVariable'><span class="jlbinding">UncertaintyQuantification.IntervalVariable</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
IntervalVariable(lb::Real, ub::Real, name::Symbol)
```


Defines an interval variable with a lower bound `lb`, upper bound `ub`, and an identifying `name`.

`IntervalVariable` can be passed directly to analyses and simulations. For other uses, such as building probability boxes (p-boxes) from interval parameters, use [`Interval`](/api/inputs#UncertaintyQuantification.Interval) instead.

**Fields**
- `lb::Real`: Lower bound of the interval.
  
- `ub::Real`: Upper bound of the interval.
  
- `name::Symbol`: Name or identifier for the variable.
  

**Examples**

```julia
julia> IntervalVariable(0.10, 0.14, :x)
x ∈ [0.1, 0.14]
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/imprecise/interval.jl#L47-L66" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.EmpiricalDistribution' href='#UncertaintyQuantification.EmpiricalDistribution'><span class="jlbinding">UncertaintyQuantification.EmpiricalDistribution</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
EmpiricalDistribution(x::Vector{<:Real}, n::Integer=10000)

Creates an empirical distribution from the data given in `x` using kernel density estimation.
The kernel used is Gaussian and the bandwidth is obtained through the Sheather-Jones method.
The support is inferred from the kde using numerical root finding.
The `cdf` and `quantile` functions are linearly interpolated using `n` data points.
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/empiricaldistribution.jl#L1-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.JointDistribution' href='#UncertaintyQuantification.JointDistribution'><span class="jlbinding">UncertaintyQuantification.JointDistribution</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
JointDistribution{D<:Union{Copula,MultivariateDistribution}, M<:Union{RandomVariable,Symbol}}(d, m)
```


Represents a joint probability distribution, either via a copula and a vector of marginal random variables, or a multivariate distribution and a vector of variable names.

**Constructors**
- JointDistribution(d::Copula, m::Vector{RandomVariable}):
  - Use a copula `d` to combine the marginal distributions in `m` into a joint distribution.
    
  - The copula&#39;s dimension must match the length of `m`.
    
  - `m` must be a vector of `RandomVariable`.
    
  
- JointDistribution(d::MultivariateDistribution, m::Vector{Symbol}):
  - Use a multivariate distribution `d` with named components specified by `m`.
    
  - The distribution&#39;s dimension (number of variables) must match the length of `m`.
    
  - `m` must be a vector of `Symbol`.
    
  

**Examples**

```julia
julia> JointDistribution(GaussianCopula([1.0 0.71; 0.71 1.0]), [RandomVariable(Normal(), :x), RandomVariable(Uniform(), :y)])
JointDistribution{Copula, RandomVariable}(GaussianCopula([1.0 0.71; 0.71 1.0]), RandomVariable[RandomVariable{Normal{Float64}}(Normal{Float64}(μ=0.0, σ=1.0), :x), RandomVariable{Uniform{Float64}}(Uniform{Float64}(a=0.0, b=1.0), :y)])
```


```julia
julia> JointDistribution(MultivariateNormal([1.0 0.71; 0.71 1.0]), [:x, :y])
JointDistribution{MultivariateDistribution, Symbol}(ZeroMeanFullNormal(
dim: 2
μ: Zeros(2)
Σ: [1.0 0.71; 0.71 1.0]
)
, [:x, :y])
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/jointdistribution.jl#L1-L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Functions {#Functions}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.sample' href='#UncertaintyQuantification.sample'><span class="jlbinding">UncertaintyQuantification.sample</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
sample(rv::RandomVariable, n::Integer=1)
```


Generates n samples from a random variable. Returns a DataFrame.

**Examples**

See also: [`RandomVariable`](/api/inputs#UncertaintyQuantification.RandomVariable)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/randomvariables/randomvariable.jl#L21-L29" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.sample-2' href='#UncertaintyQuantification.sample-2'><span class="jlbinding">UncertaintyQuantification.sample</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
sample(inputs::Vector{<:UQInput}, n::Integer=1)
```


Generates n correlated samples from a collection of inputs. Returns a DataFrame

See also: [`RandomVariable`](/api/inputs#UncertaintyQuantification.RandomVariable), [`Parameter`](/api/inputs#UncertaintyQuantification.Parameter)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/inputs.jl#L1-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.GaussianMixtureModel' href='#UncertaintyQuantification.GaussianMixtureModel'><span class="jlbinding">UncertaintyQuantification.GaussianMixtureModel</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
GaussianMixtureModel(
    data::DataFrame,
    number_components::Integer;
    maximum_iterations::Integer=100,
    tolerance::Number=1e-4,
)
```


Fits a Gaussian mixture model to the given data using the Expectation-Maximization (EM) algorithm.

**Arguments**
- `data::DataFrame`: The input data as a DataFrame, where each column represents a variable.
  
- `number_components::Integer`: The number of Gaussian components to fit.
  
- `maximum_iterations::Integer=100`: (Optional) The maximum number of EM iterations. Default is 100.
  
- `tolerance::Number=1e-4`: (Optional) The convergence tolerance for the EM algorithm. Default is 1e-4.
  

**Returns**
- `JointDistribution`: A joint distribution object representing the fitted Gaussian mixture model, with variable names corresponding to the columns of `data`.
  

See also: [`JointDistribution`](/api/inputs#UncertaintyQuantification.JointDistribution)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/gaussianmixtures.jl#L1-L22" target="_blank" rel="noreferrer">source</a></Badge>

</details>

