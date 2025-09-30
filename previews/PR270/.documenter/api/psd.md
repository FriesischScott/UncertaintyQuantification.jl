
# Power Spectral Density Functions {#Power-Spectral-Density-Functions}

Construction and evaluation of different Power Spectral Density (PSD) functions. Correpsonding theory and literature can be found here [Semi-empirical-PSD-functions](/manual/dynamics#Semi-empirical-PSD-functions).

## Index {#Index}
- [`UncertaintyQuantification.CloughPenzien`](#UncertaintyQuantification.CloughPenzien-Tuple{AbstractVector{<:Real},%20Vararg{Real,%205}})
- [`UncertaintyQuantification.EmpiricalPSD`](#UncertaintyQuantification.EmpiricalPSD)
- [`UncertaintyQuantification.KanaiTajimi`](#UncertaintyQuantification.KanaiTajimi-Tuple{AbstractVector{<:Real},%20Real,%20Real,%20Real})
- [`UncertaintyQuantification.ShinozukaDeodatis`](#UncertaintyQuantification.ShinozukaDeodatis-Tuple{AbstractVector{<:Real},%20Real,%20Real})


## Types of PSD functions {#Types-of-PSD-functions}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.CloughPenzien-Tuple{AbstractVector{<:Real}, Vararg{Real, 5}}' href='#UncertaintyQuantification.CloughPenzien-Tuple{AbstractVector{<:Real}, Vararg{Real, 5}}'><span class="jlbinding">UncertaintyQuantification.CloughPenzien</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
CloughPenzien(ω::AbstractVector{<:Real}, S_0::Real, ω_f::Real, ζ_f::Real, ω_g::Real, ζ_g::Real)
```


Constructs a `CloughPenzien` instance representing a power spectral density function with the given parameters.

**Arguments / Parameters**
- `ω::AbstractVector{<:Real}`: A vector of angular frequencies.
  
- `S_0::Real`: A scaling factor.
  
- `ω_f::Real`: Frequency parameter for the first oscillator.
  
- `ζ_f::Real`: Damping ratio for the first oscillator.
  
- `ω_g::Real`: Frequency parameter for the second oscillator.
  
- `ζ_g::Real`: Damping ratio for the second oscillator.
  

**Returns**

A discretized `CloughPenzien` power spectral density function specified by given arguments (parameters).

**Example**

```julia
w = 0:0.1:10
S_0 = 1.0
ω_f = 2.0
ζ_f = 0.05
ω_g = 3.0
ζ_g = 0.1
cp_psd = CloughPenzien(w, S_0, ω_f, ζ_f, ω_g, ζ_g)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/dynamics/psd.jl#L13-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.KanaiTajimi-Tuple{AbstractVector{<:Real}, Real, Real, Real}' href='#UncertaintyQuantification.KanaiTajimi-Tuple{AbstractVector{<:Real}, Real, Real, Real}'><span class="jlbinding">UncertaintyQuantification.KanaiTajimi</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
KanaiTajimi(ω::AbstractVector{<:Real}, S_0::Real, ω_0::Real, ζ::Real) -> KanaiTajimi
```


Constructs a `KanaiTajimi` instance representing a power spectral density function with the given parameters.

**Arguments**
- `ω::AbstractVector{<:Real}`: A vector of angular frequencies.
  
- `S_0::Real`: A scaling factor.
  
- `ω_0::Real`: Natural frequency of the oscillator.
  
- `ζ::Real`: Damping ratio of the oscillator.
  

**Returns**

A discretized `KanaiTajimi` power spectral density function specified by given arguments (parameters).

**Example**

```julia
w = 0:0.1:10
S_0 = 1.0
ω_0 = 2.0
ζ = 0.05
kt = KanaiTajimi(w, S_0, ω_0, ζ)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/dynamics/psd.jl#L60-L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.ShinozukaDeodatis-Tuple{AbstractVector{<:Real}, Real, Real}' href='#UncertaintyQuantification.ShinozukaDeodatis-Tuple{AbstractVector{<:Real}, Real, Real}'><span class="jlbinding">UncertaintyQuantification.ShinozukaDeodatis</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ShinozukaDeodatis(ω::AbstractVector{<:Real}, σ::Real, b::Real)
```


Constructs a `ShinozukaDeodatis` instance representing a power spectral density function with the given parameters.

**Arguments**
- `ω::AbstractVector{<:Real}`: A vector of angular frequencies.
  
- `σ::Real`: A hyperparamter related to the variance of the stochastic process.
  
- `b::Real`: A parameter related to the correlation length of the stochastic process.
  

**Returns**

A discretized `ShinozukaDeodatis` instance with the power spectral density function specified by given arguments (parameters).

**Example**

```julia
w = 0:0.1:10
σ = 1.0
b = 0.5
sd = ShinozukaDeodatis(w, σ, b)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/dynamics/psd.jl#L109-L129" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.EmpiricalPSD' href='#UncertaintyQuantification.EmpiricalPSD'><span class="jlbinding">UncertaintyQuantification.EmpiricalPSD</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
EmpiricalPSD(ω::AbstractVector{<:Real}, p::AbstractVector{<:Real}) -> EmpiricalPSD
```


Constructs an `EmpiricalPSD` instance with the given angular frequencies and manually provided power spectral density values.

**Arguments**
- `ω::AbstractVector{<:Real}`: A vector of angular frequencies.
  
- `p::AbstractVector{<:Real}`: A vector of power spectral density values corresponding to the frequencies in `ω`.
  

**Returns**

A discretized `EmpiricalPSD` instance with manually pre-specified provided power spectral density values.

**Example**

```julia
w = 0:0.1:10
p_values = rand(length(w))  # Example empirical PSD values
emp_psd = EmpiricalPSD(w, p_values)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/dynamics/psd.jl#L139-L157" target="_blank" rel="noreferrer">source</a></Badge>

</details>

