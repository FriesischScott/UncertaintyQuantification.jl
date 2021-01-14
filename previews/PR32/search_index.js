var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = UncertaintyQuantification\nDocTestSetup = quote\n    using UncertaintyQuantification\nend","category":"page"},{"location":"api/#Inputs","page":"API","title":"Inputs","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Parameter(value::Real, name::Symbol)\nRandomVariable(dist::Sampleable{Univariate}, name::Symbol)\nsample(inputs::Array{<:UQInput}, n::Integer)\nsample(rv::RandomVariable, n::Integer)\n","category":"page"},{"location":"api/#UncertaintyQuantification.sample-Tuple{Array{var\"#s3\",N} where N where var\"#s3\"<:UncertaintyQuantification.UQInput,Integer}","page":"API","title":"UncertaintyQuantification.sample","text":"sample(inputs::Array{<:UQInput}, n::Integer)\n\nGenerates n correlated samples from a collection of inputs. Returns a DataFrame\n\nSee also: RandomVariable, Parameter\n\n\n\n\n\n","category":"method"},{"location":"api/#UncertaintyQuantification.sample-Tuple{RandomVariable,Integer}","page":"API","title":"UncertaintyQuantification.sample","text":"sample(rv::RandomVariable, n::Integer)\n\nGenerates n samples from a random variable. Returns a DataFrame\n\nSee also: RandomVariable\n\n\n\n\n\n","category":"method"},{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Pages = [\"api.md\"]\nModule = [\"UncertaintyQuantification\"]","category":"page"},{"location":"Sensitivity/#Sensitivity-Analysis","page":"Sensitivity Analysis","title":"Sensitivity Analysis","text":"","category":"section"},{"location":"Sensitivity/","page":"Sensitivity Analysis","title":"Sensitivity Analysis","text":"TODO","category":"page"},{"location":"Sensitivity/","page":"Sensitivity Analysis","title":"Sensitivity Analysis","text":"Simulation methods\nDifferent indices","category":"page"},{"location":"UncProp/#Uncertainty-Propagation","page":"Uncertainty Propagation","title":"Uncertainty Propagation","text":"","category":"section"},{"location":"UncProp/","page":"Uncertainty Propagation","title":"Uncertainty Propagation","text":"TODO","category":"page"},{"location":"UncProp/","page":"Uncertainty Propagation","title":"Uncertainty Propagation","text":"Defining random variables\nDefining correlations\nSimulation methods\nPlotting output","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"An integrated example","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"TODO","category":"page"},{"location":"Reliability/#Reliability-Analysis","page":"Reliability Analysis","title":"Reliability Analysis","text":"","category":"section"},{"location":"Reliability/","page":"Reliability Analysis","title":"Reliability Analysis","text":"TODO","category":"page"},{"location":"Reliability/","page":"Reliability Analysis","title":"Reliability Analysis","text":"Defining failure regions\nMonte Carlo\nLine sampling","category":"page"},{"location":"#UncertaintyQuantification.jl","page":"Home","title":"UncertaintyQuantification.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia package for uncertainty quantification. Current (very limited) functionality includes:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Simulation-based Reliability Analysis (Monte Carlo, Line Sampling)\nSensitivity Analysis (local)","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ADD\nADD","category":"page"},{"location":"#Collaborators","page":"Home","title":"Collaborators","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ADD\nADD","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Two ways to install and use:","category":"page"},{"location":"","page":"Home","title":"Home","text":"1. From the julia package manager","category":"page"},{"location":"","page":"Home","title":"Home","text":"You may download the lastest release by:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]\n(v1.0) pkg> add UncertaintyQuantification\njulia> using UncertaintyQuantification","category":"page"},{"location":"","page":"Home","title":"Home","text":"or the latest version by:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]\n(v1.0) pkg> add UncertaintyQuantification#master\njulia> using UncertaintyQuantification","category":"page"},{"location":"","page":"Home","title":"Home","text":"or","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]\n(v1.0) pkg> add https://github.com/FriesischScott/UncertaintyQuantification.jl\njulia> using UncertaintyQuantification","category":"page"},{"location":"","page":"Home","title":"Home","text":"2. Downloading the source code","category":"page"},{"location":"","page":"Home","title":"Home","text":"git clone https://github.com/FriesischScott/UncertaintyQuantification.jl\n\njulia> include(\"UncertaintyQuantification.jl/src/UncertaintyQuantification.jl\")\njulia> using Main.UncertaintyQuantification","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#related-packages:","page":"Home","title":"related packages:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"OpenCossan: Matlab-based toolbox for uncertainty quantification and management\nMonteCarloMeasurements.jl: uncertainty propagation with Monte Carlo and particle-filtering\nPolyChaos.jl: polynomial chaos in Julia\nProbabilityBoundsAnalysis.jl: rigorous computations with probabilities\nIntervalArithmetic.jl: rigorous computations with intervals","category":"page"},{"location":"#Bibliography","page":"Home","title":"Bibliography","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For Uncertainty Quantification","category":"page"},{"location":"","page":"Home","title":"Home","text":"ADD\nADD","category":"page"},{"location":"","page":"Home","title":"Home","text":"For Sensitivity Analysis","category":"page"},{"location":"","page":"Home","title":"Home","text":"ADD\nADD","category":"page"},{"location":"","page":"Home","title":"Home","text":"For Reliability Analysis","category":"page"},{"location":"","page":"Home","title":"Home","text":"* De Angelis, Marco, Edoardo Patelli, and Michael Beer. \"Advanced line sampling for efficient robust reliability analysis.\" Structural safety 52 (2015): 170-182.\nADD","category":"page"},{"location":"","page":"Home","title":"Home","text":"For probability boxes:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Ferson, S., V. Kreinovich, L. Ginzburg, K. Sentz and D.S. Myers. 2003. Constructing probability boxes and Dempster-Shafer structures. Sandia National Laboratories, SAND2002-4015, Albuquerque, New Mexico\nBeer, M., S. Ferson, and V. Kreinovich. 2013. Imprecise probabilities in engineering analyses. Mechanical Systems and Signal Processing 37: 429","category":"page"},{"location":"","page":"Home","title":"Home","text":"For dependency modelling:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Ferson, S., R. Nelsen, J. Hajagos, D. Berleant, J. Zhang, W.T. Tucker, L. Ginzburg and W.L. Oberkampf. 2004. Dependence in Probabilistic Modeling, Dempster-Shafer Theory, and Probability Bounds Analysis. Sandia National Laboratories, SAND2004-3072, Albuquerque, NM\nNelsen, Roger B. An introduction to copulas. Springer Science & Business Media, 2007","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"}]
}
