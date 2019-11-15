module UncertaintyQuantification

using LinearAlgebra, DataFrames, FiniteDifferences, Dierckx, Reexport

@reexport using Distributions

import Base: rand, names

abstract type UQtypes end
abstract type AbstractInput <: UQtypes end
abstract type AbstractModel <: UQtypes end

export
    # inputs
      Parameter,
      RandomVariable,
      RandomVariableSet,

      Model,

      LineSampling,
      MonteCarlo,

    # methods
      evaluate,
      rand,
      sample,
      mean,
      gradient,
      gradient_in_standard_normal_space,
      to_standard_normal_space!,
      to_physical_space!,
      probability_of_failure

include("inputs/inputs.jl")
include("inputs/parameter.jl")
include("inputs/randomvariable.jl")
include("inputs/randomvariableset.jl")

include("models/model.jl")

include("sensitivity/gradient.jl")

include("simulations/linesampling.jl")
include("simulations/montecarlo.jl")

include("reliability/probabilityoffailure.jl")

end
