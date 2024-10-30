#===
# Bayesian Updating

## Inverse eigenvalue problem with maximum likelihood and maximum a posteriori point estimation

The inverse eigenvalue problem can also be solved with

===#

#md using UncertaintyQuantification # hide
#md using Plots # hide

#jl  using UncertaintyQuantification
#jl using Plots

Y = [
    1.51 0.33
    4.01 0.3
    3.16 0.27
    3.21 0.18
    2.19 0.33
    1.71 0.23
    2.73 0.21
    5.51 0.2
    1.95 0.11
    4.48 0.2
    1.43 0.16
    2.91 0.26
    3.81 0.23
    3.58 0.25
    2.62 0.25
]

λ1 = @. Model(df -> ((df.θ1 + 2 * df.θ2) + sqrt(df.θ1^2 + 4 * df.θ2^2)) / 2, :λ1)
λ2 = @. Model(df -> ((df.θ1 + 2 * df.θ2) - sqrt(df.θ1^2 + 4 * df.θ2^2)) / 2, :λ2)

σ = [1.0 0.1]
function likelihood(df)
    λ = [df.λ1 df.λ2]

    return log.(exp.([-0.5 * sum(((Y .- λ[n, :]') ./ σ) .^ 2) for n in axes(λ, 1)]))
end

# We will solve this problem using the TMCMC algorithm, as well as multi-objective maximum a priori (MAP) and maximum likelihood (ML) estimates. Therefore, the next step is to define the [`RandomVariable`](@ref) vector of the prior, followed by the objects for the estimaters ([`TransitionalMarkovChainMonteCarlo`](@ref), [`MaximumAPosteriori`](@ref), [`MaximumLikelihood`](@ref) ). We also have to choose number of samples and burn-in for TMCMC, as well as starting points for the optimization in ML and MAP.

prior = RandomVariable.(Normal(1, .5), [:θ1, :θ2])
burnin = 0

n = 1000

x0 = [[1., 1.],[3.,.5],[2.,2.]]

tmcmc = TransitionalMarkovChainMonteCarlo(prior, n, burnin)
MAP = MaximumAPosterioriBayesian(prior, "LBFGS", x0)
MLE = MaximumLikelihoodBayesian(prior, "LBFGS", x0)

# With the prior, likelihood, models and  MCMC sampler defined, the last step is to call the [`bayesinupdating`](@ref) method.

samples, evidence = bayesianupdating(likelihood, [λ1, λ2], tmcmc)
MapEstimate = bayesianupdating(likelihood, [λ1, λ2], MAP)
MLEstimate = bayesianupdating(likelihood, [λ1, λ2], MLE)

scatter(samples.θ1, samples.θ2; lim=[0, 4], label="TMCMC", xlabel="θ1", ylabel="θ2")
scatter!((MapEstimate.θ1, MapEstimate.θ2), label="MAP")
scatter!((MLEstimate.θ1, MLEstimate.θ2), label="MLE")

#  A scatter plot of the resulting samples shows convergence to two distinct regions. Unlike the transitional Markov Chain Monte Carlo algorithm, the standard Metropolis-Hastings algorithm would have only identified one of the two regions. Since we used a uniform prior distribution, ML and MAP estimates find the same estimates.
