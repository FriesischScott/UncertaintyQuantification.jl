using DataFrames
using Plots
using UncertaintyQuantification

n1, n2 = 200, 500
N₁ = MvNormal([2.0, 2.0], [0.5 0.0; 0.0 0.5])
N₂ = MvNormal([5.0, 3.0], [1.0 0.8; 0.8 1.5])
X = transpose([rand(N₁, n1) rand(N₂, n2)])

df = DataFrame(X, [:x1, :x2])

gmm = GaussianMixtureModel(df, 2)

x_range = range(-2, 10, length=100)
y_range = range(-2, 10, length=100)

scatter(df.x1, df.x2, alpha=0.3, label="Data")
Z = [pdf(gmm, [x, y]) for x in x_range, y in y_range]
contour!(x_range, y_range, Z', levels=10, linewidth=2, c=2, label="GMM")

samples = sample(gmm, 500)
scatter!(samples.x1, samples.x2, alpha=0.3, c=2, label="Samples")

m = MixtureModel(MvNormal[
    MvNormal([-1.0, 2.0], [1.0 0.5; 0.5 1.0]),
    MvNormal([2.0, -1.0], [1.5 0.3; 0.3 1.5]),
    MvNormal([3.0, 3.0], [1.0 0.2; 0.2 1.0])], [0.4, 0.3, 0.3]
)

gmm = JointDistribution(m, [:x1, :x2])

samples = sample(gmm, 1000)

x_range = range(-5, 8, length=100)
y_range = range(-5, 8, length=100)

scatter(samples.x1, samples.x2, alpha=0.3, label="Samples")
Z = [pdf(gmm, [x, y]) for x in x_range, y in y_range]
contour!(x_range, y_range, Z', levels=10, linewidth=2, c=2, label="GMM")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
