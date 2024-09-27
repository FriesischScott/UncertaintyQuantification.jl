using DataFrames
using UncertaintyQuantification
using AbstractGPs

l = Parameter(1.8, :l) # length
b = Parameter(0.12, :b) # width

h = RandomVariable(Normal(0.24, 0.01), :h) # height

μ, σ = distribution_parameters(10e9, 1.6e9, LogNormal)
E = RandomVariable(LogNormal(μ, σ), :E) # young's modulus

μ, σ = distribution_parameters(5000, 400, LogNormal)
P = RandomVariable(LogNormal(μ, σ), :P) # tip load

μ, σ = distribution_parameters(600, 140, LogNormal)
ρ = RandomVariable(LogNormal(μ, σ), :ρ) # density

c = GaussianCopula([1 0.8; 0.8 1])
jd = JointDistribution([E, ρ], c)

inputs = [l, b, h, P, jd]

inertia = Model(df -> df.b .* df.h .^ 3 / 12, :I)

displacement = Model(
    df ->
        (df.ρ .* 9.81 .* df.b .* df.h .* df.l .^ 4) ./ (8 .* df.E .* df.I) .+
        (df.P .* df.l .^ 3) ./ (3 .* df.E .* df.I),
    :w,
)

df = UncertaintyQuantification.sample(inputs, 20)
evaluate!([inertia, displacement], df)

random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
random_names = names(random_inputs)
output = :w

IN = InputNormalizer(df, random_names, true)
Y = Matrix(df[:, random_names])
Y_N = IN(df)

UQIN = UQInputNormalizer(random_inputs, true)
UQY = Matrix(df[:, random_names])
UQY_N = UQIN(df)

# Generate toy data.
num_dims_in = 3
num_dims_out = 2
num_obs = 100
X = randn(num_obs, num_dims_in)
Y = randn(num_obs, num_dims_out)

# Convert to format required for AbstractGPs / KernelFunctions.
# See docstrings for more info. This is basically a no-op.
x, y = prepare_isotopic_multi_output_data(RowVecs(X), RowVecs(Y))

# Construct multi-output model.
f = GP(LinearMixingModelKernel([SEKernel(), Matern52Kernel()], randn(2, num_dims_out)))

# Do the usual things that you would do with a single-output GP.
fx = f(x, 0.5)
logpdf(fx, y)
y_from_prior = rand(fx)
fx_mean = mean(fx)

f_post = posterior(fx, Y[:])

mean(f_post(x))
f_post(x)

y_shaped = reshape(y_from_prior, :, num_dims_out)