module UncertaintyQuantification

using Bootstrap
using DataFrames
using Dates
using Dierckx
using Distributed
using DynamicPolynomials
using FastGaussQuadrature
using FiniteDifferences
using Formatting
using GaussianProcesses
using KernelDensity
using LinearAlgebra
using Mustache
using Primes
using QuasiMonteCarlo
using Random
using Reexport
using StatsBase

@reexport using Distributions

import Base: rand, names, copy, run, length
import Distributions: cdf, quantile, pdf, logpdf, minimum, maximum, insupport
import Statistics: mean, var
import Distributions: logpdf, pdf, cdf, quantile, minimum, maximum, insupport, mean, var

abstract type UQType end

abstract type UQInput <: UQType end
abstract type UQModel <: UQType end

abstract type DeterministicUQInput <: UQInput end
abstract type RandomUQInput <: UQInput end

abstract type Copula <: UQType end

abstract type AbstractMonteCarlo end
abstract type AbstractQuasiMonteCarlo <: AbstractMonteCarlo end

abstract type AbstractDesignOfExperiments end

# Types
export AbstractDesignOfExperiments
export AbstractMonteCarlo
export AbstractQuasiMonteCarlo
export Copula
export DeterministicUQInput
export RandomUQInput
export UQInput
export UQModel
export UQType

# Structs
export EmpiricalDistribution
export BackwardFiniteDifferences
export BoxBehnken
export CentralComposite
export CentralFiniteDifferences
export ExternalModel
export Extractor
export FaureSampling
export FORM
export ForwardFiniteDifferences
export FractionalFactorial
export FullFactorial
export GaussianCopula
export GaussQuadrature
export HaltonSampling
export HermiteBasis
export ImportanceSampling
export JointDistribution
export LatinHypercubeSampling
export LatticeRuleSampling
export LeastSquares
export LegendreBasis
export LineSampling
export Model
export MonteCarlo
export ParallelModel
export Parameter
export PlackettBurman
export PolynomialChaosBasis
export PolynomialChaosExpansion
export PolyharmonicSpline
export RandomVariable
export ResponseSurface
export SobolSampling
export Solver
export SubSetInfinity
export SubSetInfinityAdaptive
export SubSetSimulation
export TwoLevelFactorial

# Methods
export calc
export count_rvs
export dimensions
export distribution_parameters
export doe_samples
export evaluate
export evaluate!
export gradient
export gradient_in_standard_normal_space
export mean
export multivariate_indices
export polynomialchaos
export probability_of_failure
export qmc_samples
export quadrature_nodes
export quadrature_weights
export rand
export sample
export sobolindices
export to_copula_space
export to_physical_space!
export to_standard_normal_space
export to_standard_normal_space!

include("inputs/empiricaldistribution.jl")
include("inputs/inputs.jl")
include("inputs/parameter.jl")
include("inputs/randomvariables/randomvariable.jl")
include("inputs/randomvariables/distributionparameters.jl")
include("inputs/copulas/gaussian.jl")
include("inputs/jointdistribution.jl")

include("solvers/solver.jl")
include("solvers/extractor.jl")

include("models/externalmodel.jl")
include("models/model.jl")
include("models/polyharmonicspline.jl")
include("models/responsesurface.jl")

include("models/pce/pcebases.jl")
include("models/pce/polynomialchaosexpansion.jl")

include("sensitivity/finitedifferences.jl")
include("sensitivity/gradient.jl")

include("simulations/doe.jl")
include("simulations/linesampling.jl")
include("simulations/montecarlo.jl")
include("simulations/subset.jl")

include("reliability/form.jl")
include("simulations/importancesampling.jl")
include("reliability/probabilityoffailure.jl")
include("sensitivity/sobolindices.jl")

include("util/wrap.jl")

end
