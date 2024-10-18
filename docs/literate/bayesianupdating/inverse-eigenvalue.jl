#===
# Bayesian Updating

## Inverse eigenvalue problem

The inverse eigenvalue problem is a classic engineering example. Here we will use Bayesian updating to sample from a bivariate posterior distribution describing unknown quantities of a matrix

```math
\begin{bmatrix}
 \theta_1 + \theta_2 & -\theta_2 \\
 -\theta_2 & \theta_2 \\
\end{bmatrix}
```

A matrix of this form can represent different problems, like the stiffness matrix describing a tuned mass damper system. In this example we assume the fixed values $\theta_1$ = 0.5 and $\theta_2$ = 1.5 for the variables.

The eigenvalues $\lambda_1$ and $\lambda_2$ of this matrix represent a physical measurable property corrupted by "noise" created for example due to environmental factors or measurement inaccuracy.

```math
\lambda_1^{noisy} = \frac{(\theta_1+2\theta_2)+\sqrt{\theta_1^2+4{\theta_2}^2}}{2} + \epsilon_1
```

```math
\lambda_2^{noisy} = \frac{(\theta_1+2\theta_2)-\sqrt{\theta_1^2+4{\theta_2}^2}}{2} + \epsilon_2
```

The "noise" terms $\epsilon_1$ and $\epsilon_2$ follow a Normal distribution with zero mean and standard deviations $\sigma_1$ = $1.0$ and $\sigma_2$ = $0.1$.

 The synthetic "noisy" data used for the Bayesian updating procedure is given in the following table.

|$\lambda_1$|$\lambda_2$|
|-----------|-----------|
|  1.51     |    0.33   |
|  4.01     |    0.30   |
|  3.16     |    0.17   |
|  3.21     |    0.18   |
|  2.19     |    0.32   |
|  1.71     |    0.23   |
|  2.73     |    0.21   |
|  5.51     |    0.20   |
|  1.95     |    0.11   |
|  4.48     |    0.20   |
|  1.43     |    0.16   |
|  2.91     |    0.26   |
|  3.91     |    0.23   |
|  3.58     |    0.25   |
|  2.62     |    0.25   |

The a priori knowledge of $\theta_1$ and $\theta_2$ is that they take values between 0.01 and 4. The likelihood function used for this problem is a bivariate Gaussian function with a covariance matrix $\begin{bmatrix} \sigma_1^2 & 0 \\ 0 & \sigma_2^2 \end{bmatrix}$, with off-diagonal terms equal to 0 and the diagonal terms corresponding to the variances of the respective noise terms.

```math
P(\lambda|\theta) \propto \exp \left[-\frac{1}{2}\sum_{i=1}^2\sum_{n=1}^{15} {\left(\frac{\lambda_{i,n}^{data}-\lambda_i^{model}}{\sigma_i}\right)}^2\right]
```

To begin the Bayesian model updating procedure we start by defining the data, the models for the eigenvalues (without the noise term) and the likelihood function.

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

prior = RandomVariable.(Uniform(0.01, 4), [:θ1, :θ2])

n = 1000
burnin = 0

x0 = [[1., 1.],[3.,.5],[2.,2.]]

tmcmc = TransitionalMarkovChainMonteCarlo(prior, n, burnin)
MAP = MaximumAPosteriori(prior, "LBFGS", x0)
MLE = MaximumLikelihood(prior, "LBFGS", x0)
GaussApprox = VIGaussianApproximation(prior, "LBFGS", x0)

# With the prior, likelihood, models and  MCMC sampler defined, the last step is to call the [`bayesianupdating`](@ref) method.

samples, evidence = bayesianupdating(likelihood, [λ1, λ2], tmcmc)
MapEstimate, MapValues = bayesianupdating(likelihood, [λ1, λ2], MAP)
MLEstimate, MLEValues = bayesianupdating(likelihood, [λ1, λ2], MLE)

elelele = bayesianupdating(likelihood,[λ1,λ2],GaussApprox)

MapEstimate = mapreduce(x->x, hcat, MapEstimate)'
MLEstimate = mapreduce(x->x, hcat, MLEstimate)'

scatter(samples.θ1, samples.θ2; lim=[0, 4], label="TMCMC", xlabel="θ1", ylabel="θ2")
scatter!((MapEstimate[:,1], MapEstimate[:,2]), label="MAP")
scatter!((MLEstimate[:,1], MLEstimate[:,2]), label="MLE")

#  A scatter plot of the resulting samples shows convergence to two distinct regions. Unlike the transitional Markov Chain Monte Carlo algorithm, the standard Metropolis-Hastings algorithm would have only identified one of the two regions. Since we used a uniform prior distribution, ML and MAP estimates find the same estimates.
