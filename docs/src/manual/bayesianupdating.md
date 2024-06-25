# Bayesian Updating

Bayesian updating is a method of statistical inference where Bayes' theorem is used to update the probability distributions of model parameters based on prior beliefs and available data.

## Bayes' Theorem

Bayes' theorem is defined as

```math
P(\theta|Y) = \frac{P(Y|\theta)P(\theta)}{P(Y)},
```

where $P(\theta)$ is the prior distribution, describing prior belief on $\theta$. $P(Y|\theta)$ is the likelihood function evaluating how the data Y supports our belief. This is a function of $\theta$ not $Y$. $P(\theta|Y)$ is called the posterior probability and expresses the probability distribution of the updated belief under data $Y$. The term $P(Y)$, often called the marginal likelihood, is the probability of the data. It can be calculated as the integral of the likelihood multiplied by the prior distribution over the sample space of $\theta$

```math
P(Y) = \int{}P(Y|\theta)P(\theta), d\theta{}.
```

This term serves as a normalizing constant for the posterior probability. However, as it can be difficult or even impossible to calculate it is often disregarded. Instead, only the product of likelihood and prior is used, as it is proportional to the posterior probability

```math
P(\theta|Y) \propto P(Y|\theta)P(\theta).
```

Based on this relationship, the posterior probability can be approximated without calculation of $P(Y)$ using a variety of sampling methods. Classic approaches such as rejection sampling can be inefficient, especially for multivariate cases, because of high rejection rates. Instead, Metropolis et al., proposed the use of Markov chains to increase efficiency [metropolisEquationStateCalculations1953](@cite).

## Markov Chain Monte Carlo

Markov chains are sequences of variables, where each variable is dependent on the last. In a discrete space $\Omega$ the series of random variables $\{X_1,X_2,\ldots,X_t\}$ is called a Marko chain if

```math
p(X_t=x_t|X_{t-1}=x_{t-1},\ldots,X_1=x_1) = p(X_t=x_t|X_{t-1}=x_{t-1})
```

A Markov chain is called ergodic or irreducible, when it is possible to reach each state from every other state with a positive probability. Markov chains that are ergodic and time-homogeneous, i.e. the probability between states doesn't depend on time, have a unique stationary distribution such that

```math
\pi(y) = \sum_{x\in\Omega}P(y|x)\pi(x).
```

The goal of Markov chain Monte Carlo (MCMC) sampling methods is to construct a Markov chain, whose stationary distribution is equal to the posterior distribution of Bayes' theorem. This will result in samples generated from the Markov chain being equivalent to random samples of the desired distribution. The very first MCMC algorithm is the Metropolis-Hastings (MH) Algorithm.

### Metropolis Hastings

The Metropolis-Hastings algorithm, in its generalized form, was published in 1970 by W. K. Hastings [hastingsMonteCarloSampling1970](@cite). The MH algorithm is a random-walk algorithm that provides a selection criteria for choosing the next sample (\theta_{i+1}) in a Markov chain. This is done through a so-called proposal distribution $q(\theta_{i+1}|\theta_i)$ which is well known and relatively easy to sample from. Usually, symmetric proposal distributions centred at $\theta_i$ are used which makes the Normal and Uniform distributions ideal candidates. A candidate sample $\theta^*$ is sampled from the proposal distribution and accepted with probability $\alpha$

```math
\alpha = \min\left[1,\frac{P(\theta^*|Y)}{P(\theta_i|Y)}\cdot{}\frac{q(\theta_i|\theta^*)}{q(\theta^*|\theta_i)}\right].
```

Substituting the posterior with Bayes' theorem yields

```math
\alpha = \min\left[1,\frac{P(Y|\theta^*)\cdot{}P(\theta^*)/P(Y)}{P(Y|\theta_i)\cdot{}P(\theta_i)/P(Y)}\cdot{}\frac{q(\theta_i|\theta^*)}{q(\theta^*|\theta_i)}\right].
```

Note, how the normalization constant $P(Y)$ cancels out. Because of the symmetry $q(\theta_i|\theta^*) = q(\theta^*|\theta_i)$ the acceptance probability simplifies to

```math
\alpha = \min\left[1,\frac{P(\theta^*|Y)}{P(\theta_i|Y)}\right].
```

In practice, a random number $r \sim U(0,1)$ is sampled, and the candidate sample is accepted if $a \leq r$.

As an example consider a data sequence `Y` as the outcome of 100 Bernoulli trials with unknown success probability `p` (here p=0.8).

```@example metropolis
using UncertaintyQuantification # hide
 n = 100
 Y = rand(n) .<= 0.8
 return nothing # hide
```

The likelihood function which, similar to a `Model` must accept a `DataFrame`, follows a Binomial distribution. And the prior is chosen as a beta distribution with $\alpha=\beta=1$. It is often beneficial to use the log-likelihood and log-prior for numerical reasons.

```@example metropolis
    function loglikelihood(df)
            return [
                sum(logpdf.(Binomial.(n, df_i.p), sum(Y))) for df_i in eachrow(df)
            ]
        end

logprior = df -> logpdf.(Beta(1,1), df.p)
return nothing # hide
```

**UncertaintyQuantification.jl** implements a variant of the MH algorithm known as single-component Metropolis-Hastings, where the proposal and acceptance step is performed independently for each dimension. To run the algorithm, we must first define the `SingleComponentMetropolisHastings` object which requires the `UnivariateDistribution` `proposal`, a `NamedTuple` `x0` which defines the starting point of the Markov chain, the number of samples and the number of burn-in samples. The burn-in samples are used to start the chain but later discarded.

```@example metropolis
    proposal = Normal(0, 0.2)
    x0 = (;p=0.5)
    n_samples= 2000
    burnin = 500

    mh = SingleComponentMetropolisHastings(proposal, x0, n_samples, burnin)
```

The final optional argument `islog=true` can be omitted when passing the log-likelihood and log-prior. When set to `false`, the algorithm will automatically compute the `log` for both functions. Finally, the algorithm is executed using the `bayesianupdating` function. This function returns the samples and the average acceptance rate.

```@example metropolis

mh_samples, α   = bayesianupdating(logprior, loglikelihood, mh)
```

```@example metropolis
using Plots # hide
histogram(mh_samples.p, bin = 100, label = "MH"; xlabel="p")
```

### Transitional Markov Chain Monte Carlo
