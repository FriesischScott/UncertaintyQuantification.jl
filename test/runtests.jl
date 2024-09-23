using DataFrames
using Distributed
using Formatting
using HCubature
using HypothesisTests
using InteractiveUtils
using QuasiMonteCarlo
using Random
using StatsBase: fit, Histogram, corkendall
using Test
using UncertaintyQuantification

include("inputs/empericaldistribution.jl")
include("inputs/parameter.jl")
include("inputs/jointdistribution.jl")
include("inputs/imprecise/interval.jl")
include("inputs/imprecise/p-box.jl")
include("inputs/randomvariables/randomvariable.jl")
include("inputs/randomvariables/distributionparameters.jl")
include("inputs/jointdistribution.jl");
include("inputs/inputs.jl")
include("inputs/copulas/gaussian.jl")

include("models/external/solvers.jl")
include("models/external/externalmodel.jl")
include("models/model.jl")
include("models/polyharmonicspline.jl")
include("models/pce/pcebases.jl")
include("models/pce/polynomialchaosexpansion.jl")
include("models/responsesurface.jl")
include("models/imprecise/propagation.jl")

include("modelupdating/bayesianupdating.jl")

include("reliability/form.jl")
include("reliability/probabilityoffailure.jl")
include("reliability/probabilityoffailure_imprecise.jl")

include("sensitivity/gradient.jl")
include("sensitivity/sobolindices.jl")

include("simulations/doe.jl")
include("simulations/montecarlo.jl")
include("simulations/subset.jl")

if Sys.islinux()
    include("hpc/slurm.jl")
end
