# UncertaintyQuantification.jl

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
 * Third-party solvers
   * Connect to any solver by injecting random samples into source files

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

## Authors

 * Jasper Behrensdorf, Institute for Risk and Reliability, Leibniz University Hannover
 * Ander Gray, Institute for Risk and Uncertainty, University of Liverpool


### Related packages:
* [OpenCossan](https://github.com/cossan-working-group/OpenCossan): Matlab-based toolbox for uncertainty quantification and management