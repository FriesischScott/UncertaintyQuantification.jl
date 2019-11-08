module UncertaintyQuantification

using LinearAlgebra, DataFrames, FiniteDifferences, Reexport

@reexport using Distributions

import Base: rand

abstract type UQtypes end
abstract type AbstractInput <: UQtypes end
abstract type AbstractModel <: UQtypes end

export
    # inputs
      Parameter,
      RandomVariable,
      RandomVariableSet,

      Model,

      MonteCarlo,

    # methods
      evaluate,
      rand,
      sample,
      mean,
      var,
      gradient,

      probability_of_failure

include("inputs/parameter.jl")
include("inputs/randomvariable.jl")
include("inputs/randomvariableset.jl")
include("inputs/sample.jl")

include("models/model.jl")

include("sensitivity/gradient.jl")

include("simulations/montecarlo.jl")

include("reliability/probabilityoffailure.jl")

end
