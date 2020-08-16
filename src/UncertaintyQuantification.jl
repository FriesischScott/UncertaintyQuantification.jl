module UncertaintyQuantification

using LinearAlgebra, DataFrames, FiniteDifferences, Dierckx, Reexport

@reexport using Distributions

import Base: rand, names
import Statistics: mean

abstract type UQType end

abstract type UQInput <: UQType end
abstract type UQModel <: UQType end

abstract type DeterministicUQInput <: UQInput end
abstract type RandomUQInput <: UQInput end

abstract type Copula <: UQType end

export
    # inputs
      Parameter,
      RandomVariable,
      JointDistribution,
      GaussianCopula,

      Model,
      PolyharmonicSpline,

      LineSampling,
      MonteCarlo,

    # methods
      evaluate!,
      rand,
      sample,
      count_rvs,
      dimensions,
      mean,
      gradient,
      gradient_in_standard_normal_space,
      to_standard_normal_space,
      to_standard_normal_space!,
      to_physical_space!,
      to_copula_space,
      probability_of_failure,
      sobolindices,
      calcpolyspline

include("inputs/inputs.jl")
include("inputs/parameter.jl")
include("inputs/randomvariable.jl")
include("inputs/jointdistribution.jl")

include("inputs/copulas/gaussian.jl")

include("models/model.jl")
include("models/polyharmonicspline.jl")

include("sensitivity/gradient.jl")

include("simulations/linesampling.jl")
include("simulations/montecarlo.jl")

include("reliability/probabilityoffailure.jl")
include("sensitivity/sobolindices.jl")

end
