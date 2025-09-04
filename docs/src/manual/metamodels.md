# Metamodels

## Design Of Experiments

Design Of Experiments (DOE) offers various designs that can be used for creating a model of a given system. The core idea is to evaluate significant points of the system in order to obtain a sufficient model while keeping the effort to achieve this relatively low. Depending on the parameters, their individual importance and interconnections, different designs may be adequate.

The ones implemented here are `TwoLevelFactorial`, `FullFactorial`, `FractionalFactorial`, `CentralComposite`, `BoxBehnken` and `PlackettBurman`.

## Response Surface

A Response Surface is a simple polynomial surrogate model. It can be trained by providing it with evaluated points of a function or any of the aforementioned experimental designs.

## Gaussian Process Regression

### Theoretical Background
A Gaussian Process (GP) is a collection of random variables, any finite subset of which has a joint Gaussian distribution. It is fully specified by a mean function ``m(x)`` and a covariance (kernel) function ``k(x, x')``. In GP regression, we aim to model an unknown function ``f(x)``. Before observing any data, we assume that the function ``f(x)`` is distributed according to a GP:

```math
f(x) \sim \mathcal{G}\mathcal{P}\left( m(x), k(x, x')  \right).
```

This prior GP specifies that any finite collection of function values follows a multivariate normal distribution. 

#### Posterior Gaussian Process
The posterior Gaussian Process represents the distribution of functions after incorporating observed data. We denote the observation data as: 

```math
\mathcal{D} = \lbrace (\hat{x}_i, \hat{f}_i) \mid i=1, \dots, N \rbrace,
```

where ``\hat{f}_i = f(\hat{x}_i)`` in the noise-free observation case, and ``\hat{f}_i = f(\hat{x}_i) + \varepsilon_i`` in the noisy case, with independent noise terms ``\varepsilon_i \sim \mathcal{N}(0, \sigma_\varepsilon^2)``. Let ``\hat{X} = [\hat{x}_1, \dots, \hat{x}_N]`` denote the collection of observation data locations. The corresponding mean vector and covariance matrix are:

```math
\mu(\hat{X}) = [m(\hat{x}_1), \dots, m(\hat{x}_N)], \quad K(\hat{X}, \hat{X}) \text{ with entries } K_{ij} = k(\hat{x}_i, \hat{x}_j).
 ```

For a new input location ``x^*`` we are interested at the unknown function value ``f^* = f(x^*)``. By the definition of a GP, the joint distribution of observed outputs ``\hat{f}_i`` and the unknown ``f^*`` is multivariate Gaussian:

```math
\begin{bmatrix} \hat{f}\\ f^* \end{bmatrix} = \mathcal{N}\left( \begin{bmatrix} \mu(\hat{X}) \\ m(x^*) \end{bmatrix},  \begin{bmatrix} K(\hat{X}, \hat{X}) & K(\hat{X}, x^*)\\ K(x^*, \hat{X}) & K(x^*, x^*) \end{bmatrix} \right),
```

where:
- ``K(\hat{X}, \hat{X})`` is the covariance matrix with entries ``K_{ij} = k(\hat{x}_i, \hat{x}_j)``,
- ``K(\hat{X}, x^*)`` is the covariance matrix with entries ``K_{i1} = k(\hat{x}_i, x^*)``,
- and ``K(x^*, x^*)`` is the variance at the unknown input location.

We can then obtain the posterior distribution of ``f^*`` from the properties of multivariate Gaussian distributions (see, e.g. Appendix A.2 in [rasmussen2005gaussian](@cite)), by conditioning the joint Gaussian on the observed outputs ``\hat{f}_i``:

```math
f^* \mid \hat{X}, \hat{f}, x^* \sim \mathcal{N}(\mu^*(x^*), \Sigma^*(x^*)),
```

with 

```math
\mu^*(x^*) = m(x^*) + K(x^*, \hat{X})K(\hat{X}, \hat{X})^{-1}(\hat{f} - \mu(\hat{X})), \\
\Sigma^*(x^*) = K(x^*, x^*) - K(x^*, \hat{X})K(\hat{X}, \hat{X})^{-1}K(\hat{X}, x^*).
```

In the noisy observation case, the covariance between training points is adjusted by adding the noise variance::

```math
K(\hat{X}, \hat{X}) \rightarrow K(\hat{X}, \hat{X}) + \sigma^2_{\varepsilon}I.
```

The computation of the posterior predictive distribution generalizes straightforwardly to multiple input locations, providing both the posterior mean, which can serve as a regression estimate of the unknown function, and the posterior variances, which quantify the uncertainty at each point. Because the posterior is multivariate Gaussian, one can also sample function realizations at specified locations to visualize possible functions consistent with the observed data.

#### Hyperparameter optimization
The GP prior, together with the observed data, defines a posterior distribution over functions that captures predictions at new inputs, including uncertainty. The noise-free and noisy cases differ only in the posterior covariance, which incorporates the observation noise when present.

### Constructing A Gaussian Process Regression Model