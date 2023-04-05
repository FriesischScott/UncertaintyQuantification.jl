UncertaintyQuantification.jl
============================

![Build Status](https://github.com/friesischscott/UncertaintyQuantification.jl/workflows/CI/badge.svg)
[![Coverage Status](https://codecov.io/gh/FriesischScott/UncertaintyQuantification.jl/branch/master/graph/badge.svg?token=LfslMAoWvA)](https://codecov.io/gh/FriesischScott/UncertaintyQuantification.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3993816.svg)](https://doi.org/10.5281/zenodo.3993816)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://friesischscott.github.io/UncertaintyQuantification.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://friesischscott.github.io/UncertaintyQuantification.jl/dev)

A Julia package for uncertainty quantification. Current functionality includes:

* Simulation-based reliability analysis
  * Monte Carlo simulation
  * Quasi Monte Carlo simulation (Sobol, Halton, Latin Hypercube, Lattice Rule)
  * Line Sampling
  * Subset Simulation
* Sensitivity analysis
  * Gradients
  * Sobol indices
* Metamodeling
  * Polyharmonic splines
  * Polynomial Chaos Expansion
  * Response Surface Methodology
* Third-party solvers
  * Connect to any solver by injecting random samples into source files
