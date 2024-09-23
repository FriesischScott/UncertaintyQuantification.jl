# Bayesian Updating

Bayesian updating is a method of statistical inference where Bayes' theorem is used to update the probability distributions of model parameters based on prior beliefs and available data.

## Bayes' Theorem

Bayes' theorem is defined as

```math
P(\theta|Y) = \frac{P(Y|\theta)P(\theta)}{P(Y)},
```

where $P(\theta)$ is the prior distribution, describing prior belief on $\theta$. $P(Y|\theta)$ is the likelihood function evaluating how the data Y supports our belief. This is a function of $\theta$ not $Y$. $P(\theta|Y)$ is called the posterior distribution and expresses an updated belief under data $Y$. The term $P(Y)$, often called the marginal likelihood, or the evidence. It can be calculated as the integral of the likelihood multiplied by the prior distribution over the sample space of $\theta$

```math
P(Y) = \int{}P(Y|\theta)P(\theta), d\theta{}.
```

This term serves as a normalizing constant for the posterior probability. However, as it can be difficult or even impractical to calculate it is often disregarded. Instead, only the product of likelihood and prior is used, as it is proportional to the posterior probability

```math
P(\theta|Y) \propto P(Y|\theta)P(\theta).
```

Based on this relationship, the posterior probability can be approximated without calculation of $P(Y)$ using a variety of sampling methods. Classic approaches such as rejection sampling can be inefficient, especially for multivariate cases due to high rejection rates. Instead, Metropolis et al., proposed the use of Markov chains to increase efficiency [metropolisEquationStateCalculations1953](@cite).

## Markov Chain Monte Carlo

Markov chains are sequences of variables, where each variable is dependent on the last. In a discrete space $\Omega$ the series of random variables $\{X_1,X_2,\ldots,X_t\}$ is called a Marko chain if

```math
p(X_t=x_t|X_{t-1}=x_{t-1},\ldots,X_1=x_1) = p(X_t=x_t|X_{t-1}=x_{t-1}) .
```

A Markov chain is called ergodic or irreducible when it is possible to reach each state from every other state with a positive probability. Markov chains that are ergodic and time-homogeneous, i.e. the probability between states doesn't depend on time, and have a unique stationary distribution such that

```math
\pi(y) = \sum_{x\in\Omega}P(y|x)\pi(x).
```

The goal of Markov chain Monte Carlo (MCMC) sampling methods is to construct a Markov chain, whose stationary distribution is equal to the posterior distribution of Bayes' theorem. This will result in samples generated from the Markov chain being equivalent to random samples of the desired distribution. The very first MCMC algorithm is the Metropolis-Hastings (MH) Algorithm.

### Metropolis Hastings

The Metropolis-Hastings algorithm, was published in 1970 by W. K. Hastings [hastingsMonteCarloSampling1970](@cite). The MH algorithm is a random-walk algorithm that provides a selection criteria for choosing the next sample $(\theta_{i + 1})$ in a Markov chain. This is done through a so-called proposal distribution $q(\theta_{i + 1}|\theta_i)$ which is well known and relatively easy to sample from. Usually, symmetric proposal distributions centred at $(\theta_i)$ are used, for example Normal and Uniform distributions. A candidate sample $\theta^*$ is sampled from the proposal distribution and accepted with probability $\alpha$

```math
\alpha = \min\left[1,\frac{P(\theta^*|Y)}{P(\theta_i|Y)}\cdot{}\frac{q(\theta_i|\theta^*)}{q(\theta^*|\theta_i)}\right].
```

Substituting the posterior with Bayes' theorem yields

```math
\alpha = \min\left[1,\frac{P(Y|\theta^*)\cdot{}P(\theta^*)/P(Y)}{P(Y|\theta_i)\cdot{}P(\theta_i)/P(Y)}\cdot{}\frac{q(\theta_i|\theta^*)}{q(\theta^*|\theta_i)}\right].
```

Note, how the normalization constant $P(Y)$ cancels out. When the proposal is symmetric $q(\theta_i|\theta^*) = q(\theta^*|\theta_i)$ the acceptance probability further simplifies to

```math
\alpha = \min\left[1,\frac{P(\theta^*|Y)}{P(\theta_i|Y)}\right].
```

In practice, a random number $r \sim U(0,1)$ is sampled, and the candidate sample is accepted if $a \leq r$

```math
\theta_{i + 1} = \theta^*     \qquad  \text{if} \quad a \leq r,\\
\theta_{i + 1} = \theta_{i}   \qquad  \text{otherwise.} \\
```

As an example consider a synthetic data sequence `Y` as the outcome of 100 Bernoulli trials with unknown success probability `p` (here p=0.8).

```@example metropolis
using UncertaintyQuantification # hide
 n = 100
 Y = rand(n) .<= 0.8
 return nothing # hide
```

The likelihood function which, similar to a `Model` must accept a `DataFrame`, follows a Binomial distribution and returns the likelihood for each row in the `DataFrame` as a vector. The prior is chosen as a beta distribution with $\alpha=\beta=1$ (uniform on $[0, 1]$). It is often beneficial to use the log-likelihood and log-prior for numerical reasons.

```@example metropolis
    function loglikelihood(df)
            return [
                sum(logpdf.(Binomial.(n, df_i.p), sum(Y))) for df_i in eachrow(df)
            ]
        end

logprior = df -> logpdf.(Beta(1,1), df.p)
return nothing # hide
```

**UncertaintyQuantification.jl** implements a variant of the MH algorithm known as single-component Metropolis-Hastings, where the proposal and acceptance step is performed independently for each dimension. To run the algorithm, we must first define the `SingleComponentMetropolisHastings` object which requires the `UnivariateDistribution` as a `proposal`, a `NamedTuple` for `x0` which defines the starting point of the Markov chain, the number of samples and the number of burn-in samples. The burn-in samples are used to start the chain but later discarded.

```@example metropolis
    proposal = Normal(0, 0.2)
    x0 = (;p=0.5)
    n_samples= 4000
    burnin = 500

    mh = SingleComponentMetropolisHastings(proposal, x0, n_samples, burnin)
```

The final optional argument `islog=true` can be omitted when passing the log-likelihood and log-prior. When set to `false`, the algorithm will `log` for both the likelihood and prior. Finally, the algorithm is executed using the `bayesianupdating` function. This function returns the samples and the average acceptance rate.

```@example metropolis

mh_samples, α   = bayesianupdating(logprior, loglikelihood, mh)
return nothing # hide
```

The following figure shows a histogram of the samples returned by the Metropolis-Hastings algorithm. For comparison, we also plot the analytical posterior distribution obtained using [conjugate priors](https://en.wikipedia.org/wiki/Conjugate_prior) [raiffaAppliedStatisticalDecision1961](@cite).

```@example metropolis

using Plots # hide
p = histogram(mh_samples.p, normalize=:pdf, bin = 100, label = "MH"; xlabel="p") # hide

posterior = Beta(1+ sum(Y), 1 + n - sum(Y)) # hide

x = range(0,1; length=200) # hide
y = pdf.(posterior, x) # hide

plot!(x, y, label="analytical posterior", linewidth=2) # hide

xlims!(0.5, 1.0) # hide
return p # hide
```

As a second example we will attempt to sample from a bimodial target distribution in two dimensions. The prior is uniform over $[-2, 2]$ in each dimension and the likelihood is a mixture of two Gaussian functions centred at $[0.5, 0.5]$ and $[-0.5, -0.5]$.  The standard deviation for both Gaussians are identical and if small enough will effectively disconnect the two functions.

```@example tmcmc
using UncertaintyQuantification # hide
using Plots # hide
using DataFrames # hide
prior = Uniform(-2, 2)
logprior = df -> logpdf.(prior, df.x) .+ logpdf.(prior, df.y)

N1 = MvNormal([-0.5, -0.5], 0.1)

N2 = MvNormal([0.5, 0.5], 0.1)

loglikelihood =
    df -> log.([0.5 * pdf(N1, collect(x)) + 0.5 * pdf(N2, collect(x)) for x in eachrow(df)])

n = 2000
burnin = 500

x0 = (; x=0.0, y=0.0)

proposal = Normal()

mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

mh_samples, α = bayesianupdating(logprior, loglikelihood, mh)

scatter(mh_samples.x, mh_samples.y; lim=[-2, 2], legend = :none)

xs = range(-2, 2, length = 100); ys = range(-2, 2, length = 100) # hide
sample_points = reduce(vcat,[[x y] for x in xs, y in ys][:]) # hide
df_points = DataFrame(sample_points, :auto) # hide
likelihood_eval = exp.(loglikelihood(df_points)) # hide
contour!(xs, ys, likelihood_eval, lim = [-2,2], legend = :none) # hide

```

The scatter plot clearly shows that the MH algorithm has converged to only one of the two peaks of the bimodial target (contour also plotted). In fact, this is a known weakness of the MH algorithm. However, there are a number of alternative MCMC methods that aim to solve this problem. One of these methods, known as Transitional Markov Chain Monte Carlo [chingTransitionalMarkovChain2007](@cite), will be presented next.

### Transitional Markov Chain Monte Carlo

The Transitional Markov Chain Monte Carlo (TMCMC) method [chingTransitionalMarkovChain2007](@cite) is an extension of the MH algorithm where instead of directly sampling a complex posterior distribution, samples are obtained from a series of simpler **transitional** distributions. The samples are obtained from independent single-step Markov Chains. The transitional distributions are defined as

```math
    P^j \propto P(Y|\theta)^{\beta_j} \cdot P(\theta),
```

where $j \in \{1, \ldots, m \}$ is the number of the transition step and $\beta_j$ is a tempering parameter with $\beta_1 < \cdots, \beta_m =1$. This enables a slow transition from the prior to the posterior distribution. An important part of the TMCMC algorithm is that the tempering parameter has to be selected as to ensure the transition is smooth and gradual. The algorithm's authors suggest choosing the parameter such that a coefficient of variation of 100% is maintained in the likelihood $P(Y\theta_i)^{\beta_j-\beta_{j - 1}}$. At each level $j$ the starting points for the independent Markov Chains are randomly samples (with replacement) from the current set of samples using statistical weights

```math
w(\theta_i) = \frac{P(Y|\theta_i)^{\beta_j-\beta_{j-1}}}{\sum_{i=1}^N P(Y|\theta_i)^{\beta_j-\beta_{j-1}}}.
```

The complete TMCMC algorithm can be summarized as

1. Set $j=0$ and $\beta_j=0$. Sample $\theta_i \sim P(\theta)$.
2. Set $j = j+1$.
3. Compute the next tempering parameter $\beta_j$.
4. Determine the weights $w(\theta_i)$.
5. Generate a single-step Markov chain for each $\theta_i$.
6. Repeat steps (2) to (5) until (and including) $(\beta_j=1)$.

Returning to the bimodial example, this time using the TMCMC algorithm. In order to apply a different MCMC algorithm we only need to construct a `TransitionalMarkovChainMonteCarlo` object and pass it to the `bayesianupdating` method. The definition of prior and likelihood remains the same. In difference to the `SingleComponentMetropolisHastings` the log evidence is returned instead of the acceptance rate.

```@example tmcmc

tmcmc = TransitionalMarkovChainMonteCarlo(RandomVariable.(Uniform(-2,2), [:x, :y]), n, burnin)

tmcmc_samples, S = bayesianupdating(logprior, loglikelihood, tmcmc)

scatter(tmcmc_samples.x, tmcmc_samples.y; lim=[-2, 2], label="TMCMC")
```

The resulting scatter plot shows how TMCMC is able to sample both peaks of the bimodal target distribution. The standard implementation of TMCMC uses a multivariate Gaussian proposal distribution centred at each $\theta_i$ with covariance matrix $\Sigma$ estimated from the current likelihood scaled by a factor $\beta^2$. This scaling factor defaults to $0.2$ as suggested by the authors, but can optionally be passed to the constructor as a fourth argument. Application of different MCMC Algorithms nested in the TMCMC give rise to variants of the algorithm. For example, it is possible to use the previously introduced `SingleComponentMetropolisHastings` resulting in `SingleComponentTransitionalMarkovChainMonteCarlo`.

!!! note "Note"
    `SingleComponentTransitionalMarkovChainMonteCarlo` is currently not available but planned for implementation.

For convenience, the prior can be automatically constructed from the random variables passed to `TransitionalMarkovChainMonteCarlo`.

```julia
tmcmc = TransitionalMarkovChainMonteCarlo(RandomVariable.(Uniform(-2,2), [:x, :y]), n, burnin)

tmcmc_samples, S = bayesianupdating(loglikelihood, tmcmc)
```

## Bayesian calibration of computer simulations

*UncertaintyQuantification.jl* allows for complex computer models to be included in the Bayesian updating procedure, if for example one wishes to infer unknown model parameters from experimental data of model outputs. Several models can be evaluated in order to compute the likelihood function, by passing a vector of [`UQModel`](@ref)s to the [`bayesianupdating`](@ref) method. These will be executed before the likelihood is evaluated and models outputs will be available in the `DataFrame` inside the likelihood function. There is no restriction on the type of [`UQModel`](@ref) used. For example, it is possible to use an [`ExternalModel`](@ref) and call an external solver.

For a complete example refer to the [Inverse eigenvalue problem](@ref).
