module UncertaintyQuantification

using Distributions, LinearAlgebra, DataFrames

import Base: rand

#import Distributions: Normal, LogNormal


abstract type AbstractInput <:Real end
abstract type AbstractSimulation <:Real end

export
    # inputs
    Parameter,
    RandomVariable,
    RandomVariableSet,

    Model,

    MonteCarlo,

    # methods
    copularand,
    rand,
    sample,

    probabilityOfFailure,

    Normal, LogNormal

include("inputs/parameter.jl")
include("inputs/randomvariable.jl")
include("inputs/randomvariableset.jl")
include("inputs/sample.jl")

include("models/model.jl")

include("simulations/montecarlo.jl")

include("reliability/probabilityoffailure.jl")

include("util/copula.jl")

end
