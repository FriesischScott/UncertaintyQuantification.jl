#===

# Gaussian Mixture Models

Gaussian Mixture Models (GMMs) are a probabilistic model that assumes all data points are generated from a mixture of ``K`` Gaussian distributions with unknown parameters.
A GMM can be defined by its probability density function (PDF) which is a weighted sum of the individual Gaussian densities:

```math
p(x) = \sum_{k=1}^{K} \pi_k \, \mathcal{N}(x \mid \mu_k, \Sigma_k).
```

Here, ``\pi_k`` are the mixing coefficients (weights) with property ``\sum_{k=1}^K \pi_k = 1``, ``\mu_k`` are the means, and ``\Sigma_k`` are the covariance matrices of the Gaussian components.

## Fitting a Gaussian Mixture Model using the Expectation-Maximization Algorithm

One way to find the parameters of a GMM is to use the Expectation-Maximization (EM) algorithm.
Here, we show the basic steps to fit a GMM to data using the EM algorithm based on [hastieElementsStatiscialLearning2009](@cite).
The EM algorithm iteratively refines the parameters of the GMM by alternating between two steps:
1. **Expectation Step (E-step)**: Calculate the expected value of the latent variables given the current parameters.
2. **Maximization Step (M-step)**: Update the parameters to maximize the expected log-likelihood found in the E-step.

In `UncertaintyQuantification.jl`, we can construct a GMM from available data using the EM algorithm.
First we load the necessary packages to fit the GMM and visualize the results.
===#
using DataFrames
using Plots
using UncertaintyQuantification

# Then, we generate some data from two bivariate Gaussian distributions that we use to fit a GMM.
#md using Random; Random.seed!(69) # hide
n1, n2 = 200, 500
N₁ = MvNormal([2.0, 2.0], [0.5 0.0; 0.0 0.5])
N₂ = MvNormal([5.0, 3.0], [1.0 0.8; 0.8 1.5])
X = transpose([rand(N₁, n1) rand(N₂, n2)])
#md nothing # hide

# To store and process the data, we use a `DataFrame`:
df = DataFrame(X, [:x1, :x2])
#md nothing # hide

# Then, we fit a `GaussianMixtureModel` with two dimensions ($x_1$ and $x_2$) and $K=2$ components to the data stored in `df`:
gmm = GaussianMixtureModel(df, 2)

# To visually validate the fit, we can plot the data and the fitted GMM. We create a grid of points to evaluate the GMM's PDF and plot the contours.
x_range = range(-2, 10, length=100)
y_range = range(-2, 10, length=100)

scatter(df.x1, df.x2, alpha=0.3, label="Data")
Z = [pdf(gmm, [x, y]) for x in x_range, y in y_range]
contour!(x_range, y_range, Z', levels=10, linewidth=2, c=2, label="GMM")
#md savefig("data-gmm.svg"); nothing # hide

# ![Data and fitted GMM](data-gmm.svg)

# From the fitted GMM, we can also draw samples and compare them to the original data. We generate 500 samples from the GMM and plot them.
samples = sample(gmm, 500)
scatter!(samples.x1, samples.x2, alpha=0.3, c=2, label="Samples")
#md savefig("samples-gmm.svg"); nothing # hide

# ![Data and fitted GMM, and samples from the GMM](samples-gmm.svg)

#===

## Gaussian Mixture Models fitted with other packages

Alternatively, we can use other packages, such as [GaussianMixtures.jl](https://github.com/davidavdav/GaussianMixtures.jl)
to fit Gaussian Mixture Models, as the `GaussianMixtureModel` type is compatible with any `MixtureModel` type from [Distributions.jl](https://juliastats.org/Distributions.jl/stable/).

===#

# Let's start with defining a `MixtureModel` from `Distributions.jl`:
#md using Random; Random.seed!(2) # hide
m = MixtureModel(MvNormal[
    MvNormal([-1.0, 2.0], [1.0 0.5; 0.5 1.0]),
    MvNormal([2.0, -1.0], [1.5 0.3; 0.3 1.5]),
    MvNormal([3.0, 3.0], [1.0 0.2; 0.2 1.0])], [0.4, 0.3, 0.3]
)

# Then, we can create our `GaussianMixtureModel` from the `MixtureModel` and specify the names of the dimensions.
gmm = JointDistribution(m, [:x1, :x2])
#md nothing # hide

# Again, we can sample from the GMM and, for example, evaluate the PDF and plot those results:
samples = sample(gmm, 1000)

x_range = range(-5, 8, length=100)
y_range = range(-5, 8, length=100)

scatter(samples.x1, samples.x2, alpha=0.3, label="Samples")
Z = [pdf(gmm, [x, y]) for x in x_range, y in y_range]
contour!(x_range, y_range, Z', levels=10, linewidth=2, c=2, label="GMM")

#md savefig("samples-gmm-random.svg"); nothing # hide
# ![Samples from GMM using Distributions.jl](samples-gmm-random.svg)
