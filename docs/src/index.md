# UncertaintyQuantification.jl

```@raw html
---
authors:
  - name: Jane Smith
  - name: John Doe
    platform: bluesky
  - name: Lazaro
    avatar: https://avatars.githubusercontent.com/u/19525261?v=4
    platform: github
    link: https://github.com/lazarusA

---

## Authors

<Authors />

```

A Julia package for uncertainty quantification. Current functionality includes:

* Simulation-based reliability analysis
  * Monte Carlo simulation
  * Quasi Monte Carlo simulation (Sobol, Halton)
  * Line Sampling
  * Subset Simulation
* Sensitivity analysis
  * Gradients
  * Sobol indices
* Metamodeling
  * Polyharmonic splines
  * Response Surface
* Third-party solvers
  * Connect to any solver by injecting random samples into source files
  * HPC interfacing with slurm

---

## Installation

To install the latest release through the Julia package manager run:

```julia
julia> ]add UncertaintyQuantification
julia> using UncertaintyQuantification
```

or install the latest development version with:

```julia
julia> ]add UncertaintyQuantification#master
julia> using UncertaintyQuantification
```

---

### Related packages

* [OpenCossan](https://github.com/cossan-working-group/OpenCossan): Matlab-based toolbox for uncertainty quantification and management
