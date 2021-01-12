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
* [MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl): uncertainty propagation with Monte Carlo and particle-filtering
* [PolyChaos.jl](https://github.com/timueh/PolyChaos.jl): polynomial chaos in Julia
* [ProbabilityBoundsAnalysis.jl](https://github.com/JuliaIntervals/ValidatedNumerics.jl): rigorous computations with probabilities
* [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl): rigorous computations with intervals


Bibliography
---

For Uncertainty Quantification
* ADD
* ADD

For Sensitivity Analysis
* ADD
* ADD

For Reliability Analysis
* Line Sampling
* ADD

For probability boxes:
* [*Ferson, S., V. Kreinovich, L. Ginzburg, K. Sentz and D.S. Myers. 2003. Constructing probability boxes and Dempster-Shafer structures. Sandia National Laboratories, SAND2002-4015, Albuquerque, New Mexico*](https://www.osti.gov/servlets/purl/1427258)
* *Ferson, S., and J. Siegrist. 2012. Verified computation with probabilities. Pages 95-122 in Uncertainty Quantification in Scientific Computing, edited by A. Dienstfrey and R.F. Boisvert, Springer, New York*
* [*Beer, M., S. Ferson, and V. Kreinovich. 2013. Imprecise probabilities in engineering analyses. Mechanical Systems and Signal Processing 37: 429*](https://digitalcommons.utep.edu/cgi/viewcontent.cgi?article=1733&=&context=cs_techrep&=&sei-redir=1&referer=https%253A%252F%252Fscholar.google.com%252Fscholar%253Fhl%253Den%2526as_sdt%253D0%25252C5%2526q%253DBeer%25252C%252BM.%25252C%252BS.%252BFerson%25252C%252Band%252BV.%252BKreinovich.%252B2013.%252BImprecise%252Bprobabilities%252Bin%252Bengineering%252Banalyses.%252BMechanical%252BSystems%252Band%252BSignal%252BProcessing%252B37%25253A%252B429%2526btnG%253D#search=%22Beer%2C%20M.%2C%20S.%20Ferson%2C%20V.%20Kreinovich.%202013.%20Imprecise%20probabilities%20engineering%20analyses.%20Mechanical%20Systems%20Signal%20Processing%2037%3A%20429%22)

For dependency modelling:
* [*Ferson, S., R. Nelsen, J. Hajagos, D. Berleant, J. Zhang, W.T. Tucker, L. Ginzburg and W.L. Oberkampf. 2004. Dependence in Probabilistic Modeling, Dempster-Shafer Theory, and Probability Bounds Analysis. Sandia National Laboratories, SAND2004-3072, Albuquerque, NM*](https://www.osti.gov/servlets/purl/1427286)

* *Nelsen, Roger B. An introduction to copulas. Springer Science & Business Media, 2007*


---

### Acknowledgements


