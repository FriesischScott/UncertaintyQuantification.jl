module UncertaintyQuantification

using Accessors
using Bootstrap
using DataFrames
using Dates
using Dierckx
using Distributed
using FiniteDifferences
using Formatting
using LinearAlgebra
using Mustache
using QuasiMonteCarlo
using Random
using Reexport

@reexport using Distributions

import Base: rand, names, copy, run
import Statistics: mean

abstract type UQType end

abstract type UQInput <: UQType end
abstract type UQModel <: UQType end

abstract type DeterministicUQInput <: UQInput end
abstract type RandomUQInput <: UQInput end

abstract type Copula <: UQType end

abstract type AbstractMonteCarlo end
abstract type AbstractQuasiMonteCarlo <: AbstractMonteCarlo end

# Types
export AbstractMonteCarlo
export AbstractQuasiMonteCarlo
export Copula
export DeterministicUQInput
export RandomUQInput
export UQInput
export UQModel
export UQType

# Structs
export ExternalModel
export Extractor
export GaussianCopula
export HaltonSampling
export JointDistribution
export LatinHypercubeSampling
export LineSampling
export Model
export MonteCarlo
export Parameter
export PolyharmonicSpline
export RandomVariable
export SobolSampling
export Solver
export SubSetSimulation

# Methods
export calc
export count_rvs
export dimensions
export evaluate!
export gradient
export gradient_in_standard_normal_space
export mean
export probability_of_failure
export qmc_samples
export rand
export sample
export sobolindices
export to_copula_space
export to_physical_space!
export to_standard_normal_space
export to_standard_normal_space!

include("inputs/inputs.jl")
include("inputs/parameter.jl")
include("inputs/randomvariable.jl")
include("inputs/jointdistribution.jl")

include("inputs/copulas/gaussian.jl")

include("solvers/solver.jl")
include("solvers/extractor.jl")

include("models/externalmodel.jl")
include("models/model.jl")
include("models/polyharmonicspline.jl")

include("sensitivity/gradient.jl")

include("simulations/linesampling.jl")
include("simulations/montecarlo.jl")
include("simulations/subset.jl")

include("reliability/probabilityoffailure.jl")
include("sensitivity/sobolindices.jl")

end
