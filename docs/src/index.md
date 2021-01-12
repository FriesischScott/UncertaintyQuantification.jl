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
* [* De Angelis, Marco, Edoardo Patelli, and Michael Beer. "Advanced line sampling for efficient robust reliability analysis." Structural safety 52 (2015): 170-182.](https://www.sciencedirect.com/user/identity/landing?code=sprduN_D6aEaY-XtsvKU9eUVlUozUTQAiTcfKyQD&state=retryCounter%3D0%26csrfToken%3D64322e8f-2c3c-40b0-b9d5-aa74cda145eb%26idpPolicy%3Durn%253Acom%253Aelsevier%253Aidp%253Apolicy%253Aproduct%253Ainst_assoc%26returnUrl%3D%252Fscience%252Farticle%252Fpii%252FS0167473014000927%253Fcasa_token%253DS9q_6t9W2-kAAAAA%253A_z0cfyiXZs8QwNXJz788YX4LTlUCW3epsAHU8SgRH3Cr6Sd5BcreHVkWgNvp1tOiXzUMmPvQCZ4g%26uuid%3D2a45854e-3501-4e86-b720-8467969f2f3a%26prompt%3Dnone%26cid%3Darp-ad74e836-bf4a-4da1-8482-e076d81e924b)
* ADD

For probability boxes:
* [*Ferson, S., V. Kreinovich, L. Ginzburg, K. Sentz and D.S. Myers. 2003. Constructing probability boxes and Dempster-Shafer structures. Sandia National Laboratories, SAND2002-4015, Albuquerque, New Mexico*](https://www.osti.gov/servlets/purl/1427258)
* [*Beer, M., S. Ferson, and V. Kreinovich. 2013. Imprecise probabilities in engineering analyses. Mechanical Systems and Signal Processing 37: 429*](https://digitalcommons.utep.edu/cgi/viewcontent.cgi?article=1733&=&context=cs_techrep&=&sei-redir=1&referer=https%253A%252F%252Fscholar.google.com%252Fscholar%253Fhl%253Den%2526as_sdt%253D0%25252C5%2526q%253DBeer%25252C%252BM.%25252C%252BS.%252BFerson%25252C%252Band%252BV.%252BKreinovich.%252B2013.%252BImprecise%252Bprobabilities%252Bin%252Bengineering%252Banalyses.%252BMechanical%252BSystems%252Band%252BSignal%252BProcessing%252B37%25253A%252B429%2526btnG%253D#search=%22Beer%2C%20M.%2C%20S.%20Ferson%2C%20V.%20Kreinovich.%202013.%20Imprecise%20probabilities%20engineering%20analyses.%20Mechanical%20Systems%20Signal%20Processing%2037%3A%20429%22)

For dependency modelling:
* [*Ferson, S., R. Nelsen, J. Hajagos, D. Berleant, J. Zhang, W.T. Tucker, L. Ginzburg and W.L. Oberkampf. 2004. Dependence in Probabilistic Modeling, Dempster-Shafer Theory, and Probability Bounds Analysis. Sandia National Laboratories, SAND2004-3072, Albuquerque, NM*](https://www.osti.gov/servlets/purl/1427286)

* *Nelsen, Roger B. An introduction to copulas. Springer Science & Business Media, 2007*


---

### Acknowledgements


