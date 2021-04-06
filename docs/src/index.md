# UncertaintyQuantification.jl


A Julia package for uncertainty quantification. Current (very limited) functionality includes:

 * Simulation-based Reliability Analysis (Monte Carlo, Line Sampling)
 * Sensitivity Analysis (local)



### Authors

* ADD
* ADD

### Collaborators

* ADD
* ADD


---

Installation
---
Two ways to install and use:

**1. From the julia package manager**

You may download the lastest release by:
```julia
julia> ]
(v1.0) pkg> add UncertaintyQuantification
julia> using UncertaintyQuantification
```

or the latest version by:

```julia
julia> ]
(v1.0) pkg> add UncertaintyQuantification#master
julia> using UncertaintyQuantification
```

or

```julia
julia> ]
(v1.0) pkg> add https://github.com/FriesischScott/UncertaintyQuantification.jl
julia> using UncertaintyQuantification
```

**2. Downloading the source code**
```julia
git clone https://github.com/FriesischScott/UncertaintyQuantification.jl

julia> include("UncertaintyQuantification.jl/src/UncertaintyQuantification.jl")
julia> using Main.UncertaintyQuantification
```

---

### related packages:
* [OpenCossan](https://github.com/cossan-working-group/OpenCossan): Matlab-based toolbox for uncertainty quantification and management