
# Stochastic Processes (Spectral Representation) {#Stochastic-Processes-Spectral-Representation}

Stochastic process generation based on the Spectral Representation Method which utilizes Power Spectral Density Functions. Correpsonding theory and literature can be found here [Stochastic-Process-Generation](/manual/dynamics#Stochastic-Process-Generation).

## Index {#Index}
- [`UncertaintyQuantification.SpectralRepresentation`](#UncertaintyQuantification.SpectralRepresentation-Tuple{AbstractPowerSpectralDensity,%20AbstractVector{<:Real},%20Symbol})
- [`Base.names`](#Base.names-Tuple{SpectralRepresentation})
- [`UncertaintyQuantification.dimensions`](#UncertaintyQuantification.dimensions-Tuple{SpectralRepresentation})
- [`UncertaintyQuantification.evaluate`](#UncertaintyQuantification.evaluate-Tuple{SpectralRepresentation,%20AbstractVector{<:Real}})
- [`UncertaintyQuantification.sample`](#UncertaintyQuantification.sample)
- [`UncertaintyQuantification.to_physical_space!`](#UncertaintyQuantification.to_physical_space!-Tuple{SpectralRepresentation,%20DataFrame})
- [`UncertaintyQuantification.to_standard_normal_space!`](#UncertaintyQuantification.to_standard_normal_space!-Tuple{SpectralRepresentation,%20DataFrame})


## Types and Spectral Representation functions {#Types-and-Spectral-Representation-functions}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.SpectralRepresentation-Tuple{AbstractPowerSpectralDensity, AbstractVector{<:Real}, Symbol}' href='#UncertaintyQuantification.SpectralRepresentation-Tuple{AbstractPowerSpectralDensity, AbstractVector{<:Real}, Symbol}'><span class="jlbinding">UncertaintyQuantification.SpectralRepresentation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
SpectralRepresentation(psd::AbstractPowerSpectralDensity, time::AbstractVector{<:Real}, name::Symbol) -> SpectralRepresentation
```


Constructs a `SpectralRepresentation` instance representing a stochastic process generated using the spectral representation method.

**Arguments**
- `psd::AbstractPowerSpectralDensity`: An instance of a power spectral density model.
  
- `time::AbstractVector{<:Real}`: A vector of time points.
  
- `name::Symbol`: A symbol representing the name of the process.
  

**Returns**

A `SpectralRepresentation` instance with the given arguments (parameters).

**Example**

```julia
w = 0:0.1:10
S_0 = 1.0
ω_f = 2.0
ζ_f = 0.05
ω_g = 3.0
ζ_g = 0.1
psd = CloughPenzien(w, S_0, ω_f, ζ_f, ω_g, ζ_g)
t = 0:0.1:10
name = :process1
sr = SpectralRepresentation(psd, t, name)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/stochasticprocesses/spectralrepresentation.jl#L13-L40" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.sample' href='#UncertaintyQuantification.sample'><span class="jlbinding">UncertaintyQuantification.sample</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
sample(sr::SpectralRepresentation, n::Integer=1) -> DataFrame
```


Generates samples of random phase angles for a given `SpectralRepresentation` instance.

**Arguments**
- `sr::SpectralRepresentation`: An instance of the `SpectralRepresentation` struct.
  
- `n::Integer=1`: The number of samples to generate (default is 1).
  

**Returns**

A `DataFrame` containing the generated samples of random phase angles.

**Example**

```julia
w = 0:0.1:10
psd = CloughPenzien(w, 1.0, 2.0, 0.05, 3.0, 0.1)
t = 0:0.1:10
name = :process1
sr = SpectralRepresentation(psd, t, name)
samples = sample(sr, 5)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/stochasticprocesses/spectralrepresentation.jl#L69-L91" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.evaluate-Tuple{SpectralRepresentation, AbstractVector{<:Real}}' href='#UncertaintyQuantification.evaluate-Tuple{SpectralRepresentation, AbstractVector{<:Real}}'><span class="jlbinding">UncertaintyQuantification.evaluate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
evaluate(sr::SpectralRepresentation, ϕ::AbstractVector{<:Real}) -> AbstractVector{<:Real}
```


Evaluates the stochastic process for a given `SpectralRepresentation` instance and a vector of random phase angles.

**Arguments**
- `sr::SpectralRepresentation`: An instance of the `SpectralRepresentation` struct.
  
- `ϕ::AbstractVector{<:Real}`: A vector of random phase angles.
  

**Returns**

A vector of real numbers representing the evaluated stochastic process.

**Example**

```julia
w = 0:0.1:10
psd = CloughPenzien(w, 1.0, 2.0, 0.05, 3.0, 0.1)
t = 0:0.1:10
name = :process1
sr = SpectralRepresentation(psd, t, name)
ϕ = rand(Uniform(0, 2π), length(psd.ω))
process_values = evaluate(sr, ϕ)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/stochasticprocesses/spectralrepresentation.jl#L96-L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.to_standard_normal_space!-Tuple{SpectralRepresentation, DataFrame}' href='#UncertaintyQuantification.to_standard_normal_space!-Tuple{SpectralRepresentation, DataFrame}'><span class="jlbinding">UncertaintyQuantification.to_standard_normal_space!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
to_standard_normal_space!(sr::SpectralRepresentation, df::DataFrame) -> Nothing
```


Transforms the random phase angles in the given `DataFrame` from a uniform distribution to a standard normal distribution.

**Arguments**
- `sr::SpectralRepresentation`: An instance of the `SpectralRepresentation` struct.
  
- `df::DataFrame`: A `DataFrame` containing the random phase angles to be transformed.
  

**Returns**

Nothing. The `DataFrame` is modified in place.

**Example**

```julia
w = 0:0.1:10
psd = CloughPenzien(w, 1.0, 2.0, 0.05, 3.0, 0.1)
t = 0:0.1:10
name = :process1
sr = SpectralRepresentation(psd, t, name)
samples = sample(sr, 5)
to_standard_normal_space!(sr, samples)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/stochasticprocesses/spectralrepresentation.jl#L140-L163" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.to_physical_space!-Tuple{SpectralRepresentation, DataFrame}' href='#UncertaintyQuantification.to_physical_space!-Tuple{SpectralRepresentation, DataFrame}'><span class="jlbinding">UncertaintyQuantification.to_physical_space!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
to_physical_space!(sr::SpectralRepresentation, df::DataFrame) -> Nothing
```


Transforms the random phase angles in the given `DataFrame` from a standard normal distribution to a uniform distribution.

**Arguments**
- `sr::SpectralRepresentation`: An instance of the `SpectralRepresentation` struct.
  
- `df::DataFrame`: A `DataFrame` containing the random phase angles to be transformed.
  

**Returns**

Nothing. The `DataFrame` is modified in place.

**Example**

```julia
w = 0:0.1:10
psd = CloughPenzien(w, 1.0, 2.0, 0.05, 3.0, 0.1)
t = 0:0.1:10
name = :process1
sr = SpectralRepresentation(psd, t, name)
samples = sample(sr, 5)
to_standard_normal_space!(sr, samples)  # Transform to standard normal space
to_physical_space!(sr, samples)         # Transform back to physical space
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/stochasticprocesses/spectralrepresentation.jl#L171-L195" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.dimensions-Tuple{SpectralRepresentation}' href='#UncertaintyQuantification.dimensions-Tuple{SpectralRepresentation}'><span class="jlbinding">UncertaintyQuantification.dimensions</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dimensions(sr::SpectralRepresentation) -> Int
```


Returns the number of dimensions (frequencies) in the given `SpectralRepresentation` instance.

**Arguments**
- `sr::SpectralRepresentation`: An instance of the `SpectralRepresentation` struct.
  

**Returns**

An integer representing the number of dimensions (frequencies).

**Example**

```julia
w = 0:0.1:10
psd = CloughPenzien(w, 1.0, 2.0, 0.05, 3.0, 0.1)
t = 0:0.1:10
name = :process1
sr = SpectralRepresentation(psd, t, name)
num_dimensions = dimensions(sr)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/stochasticprocesses/spectralrepresentation.jl#L203-L224" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.names-Tuple{SpectralRepresentation}' href='#Base.names-Tuple{SpectralRepresentation}'><span class="jlbinding">Base.names</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
names(sr::SpectralRepresentation) -> Vector{Symbol}
```


Returns the names of the random phase angles for a given `SpectralRepresentation` instance.

**Arguments**
- `sr::SpectralRepresentation`: An instance of the `SpectralRepresentation` struct.
  

**Returns**

A vector of symbols representing the names of the random phase angles.

**Example**

```julia
w = 0:0.1:10
psd = CloughPenzien(w, 1.0, 2.0, 0.05, 3.0, 0.1)
t = 0:0.1:10
name = :process1
sr = SpectralRepresentation(psd, t, name)
phase_angle_names = names(sr)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/inputs/stochasticprocesses/spectralrepresentation.jl#L229-L250" target="_blank" rel="noreferrer">source</a></Badge>

</details>

