module UncertaintyQuantification

using Bootstrap
using CovarianceEstimation
using DataFrames
using Dates
using Dierckx
using Distributed
using DynamicPolynomials
using FastGaussQuadrature
using FiniteDifferences
using Formatting
using KernelDensity
using LinearAlgebra
using MeshAdaptiveDirectSearch
using Mustache
using Primes
using QuasiMonteCarlo
using Random
using Reexport
using StatsBase

@reexport using Distributions

import Base: rand, names, copy, run, length
import Distributions: cdf, quantile, pdf, logpdf, minimum, maximum, insupport, mean, var
import Statistics: mean, var

abstract type UQInput end
abstract type DeterministicUQInput <: UQInput end
abstract type ImpreciseUQInput <: UQInput end
abstract type RandomUQInput <: UQInput end

abstract type UQModel end

abstract type Copula end

abstract type AbstractSimulation end
abstract type AbstractMonteCarlo <: AbstractSimulation end
abstract type AbstractQuasiMonteCarlo <: AbstractMonteCarlo end

"""
    AbstractBayesianMethod

Subtypes are used to dispatch to the differenct MCMC methods in [`bayesianupdating`](@ref).

Subtypes are:

- [`SingleComponentMetropolisHastings`](@ref)
- [`TransitionalMarkovChainMonteCarlo`](@ref)
"""
abstract type AbstractBayesianMethod end
abstract type AbstractDesignOfExperiments end

abstract type AbstractHPCScheduler end

# Types
export AbstractBayesianMethod
export AbstractDesignOfExperiments
export AbstractMonteCarlo
export AbstractQuasiMonteCarlo
export AbstractSimulation
export Copula
export DeterministicUQInput
export RandomUQInput
export ImpreciseUQInput
export UQInput
export UQModel

# Structs
export AdaptiveMetropolisHastings
export EmpiricalDistribution
export BackwardFiniteDifferences
export BoxBehnken
export CentralComposite
export CentralFiniteDifferences
export DoubleLoop
export ExternalModel
export SlurmInterface
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
export Interval
export JointDistribution
export LatinHypercubeSampling
export LatticeRuleSampling
export LeastSquares
export LegendreBasis
export LineSampling
export SingleComponentMetropolisHastings
export Model
export MonteCarlo
export ParallelModel
export Parameter
export PlackettBurman
export PolynomialChaosBasis
export PolynomialChaosExpansion
export PolyharmonicSpline
export ProbabilityBox
export RandomVariable
export RandomSlicing
export ResponseSurface
export SobolSampling
export Solver
export SubSetInfinity
export SubSetInfinityAdaptive
export SubSetSimulation
export TransitionalMarkovChainMonteCarlo
export TwoLevelFactorial

# Methods
export bayesianupdating
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
export propagate_intervals!
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

include("inputs/imprecise/interval.jl")
include("inputs/imprecise/p-box.jl")

include("inputs/randomvariables/randomvariable.jl")
include("inputs/randomvariables/distributionparameters.jl")
include("inputs/copulas/gaussian.jl")
include("inputs/jointdistribution.jl")

include("solvers/solver.jl")
include("solvers/extractor.jl")

include("hpc/slurm.jl")

include("models/externalmodel.jl")
include("models/model.jl")
include("models/imprecise/propagation.jl")
include("models/polyharmonicspline.jl")
include("models/responsesurface.jl")
include("models//slicingmodel.jl")

include("models/pce/pcebases.jl")
include("models/pce/polynomialchaosexpansion.jl")

include("modelupdating/bayesianupdating.jl")

include("sensitivity/finitedifferences.jl")
include("sensitivity/gradient.jl")

include("simulations/doe.jl")
include("simulations/linesampling.jl")
include("simulations/montecarlo.jl")
include("simulations/subset.jl")

include("reliability/form.jl")
include("simulations/importancesampling.jl")
include("reliability/probabilityoffailure.jl")
include("reliability/probabilityoffailure_imprecise.jl")
include("sensitivity/sobolindices.jl")

include("util/wrap.jl")
include("util/imprecise.jl")

end
