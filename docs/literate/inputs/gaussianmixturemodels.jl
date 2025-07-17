#===

# Gaussian Mixture Models

- GMM fitted using the Expectation-Maximization (EM) algorithm, implemented based on [hastieElementsStatiscialLearning2009](@cite)
- Samples can be drawn from the fitted GMM

===#
using DataFrames
using Plots
using UncertaintyQuantification

# We start with generating some artificial data from two bivariate Gaussian distributions.
n1, n2 = 200, 500
N₁ = MvNormal([2.0, 2.0], [0.5 0.0; 0.0 0.5])
N₂ = MvNormal([5.0, 3.0], [1.0 0.8; 0.8 1.5])
X = transpose([rand(N₁, n1) rand(N₂, n2)])
#md nothing # hide

# To store and process the data, we use a `DataFrame`:
df = DataFrame(X, [:x1, :x2])
#md nothing # hide

# Then, we initialize a `GaussianMixtureModel` with two dimensions ($x_1$ and $x_2$) and $K=2$ components. Afterwards, we fit the GMM to the data using the `fit!` function.
gmm = GaussianMixtureModel([:x1, :x2], 2)
fit!(gmm, df)

# To validate the fit, we can visualize the data and the fitted GMM. We create a grid of points to evaluate the GMM's PDF and plot the contours.
x_range = range(-2, 10, length=100)
y_range = range(-2, 10, length=100)

scatter(df.x1, df.x2, alpha=0.3, label="Data")
Z = [pdf(gmm, [x, y]) for x in x_range, y in y_range]
contour!(x_range, y_range, Z', levels=10, linewidth=2, color=2, label="GMM")
#md savefig("data-gmm.svg"); nothing # hide

# ![Data and fitted GMM](data-gmm.svg)

# From the fitted GMM, we can also draw samples and compare them to the original data. We generate 500 samples from the GMM and plot them.
samples = sample(gmm, 500)
scatter!(samples.x1, samples.x2, alpha=0.5, label="Samples", c=2)
#md savefig("samples-gmm.svg"); nothing # hide

# ![Data and fitted GMM, and samples from the GMM](samples-gmm.svg)
