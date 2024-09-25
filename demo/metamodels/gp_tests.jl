using UncertaintyQuantification
using ParameterHandling
using AbstractGPs
using Random
using DataFrames

using Zygote
using Optim
using Plots

Random.seed!(20140430)

## Training data
n = 10

# For interface with random input and model
x = RandomVariable(Uniform(0, 2π), :x)
y = Model(
    df ->
        (sin.(df.x) + 0.05*randn(length(df.x))),
    :y,
)
exp_design = ExperimentalDesign(MonteCarlo(n))

# For interface with DataFrame
df = sample(x, n)
evaluate!(y, df) 

## Set up mean, kernel and noise
# mean
mean_params = (;)
mZero(θ) = ZeroMean() #Zero mean function

# kernel
kernel_params = (;
    σ = positive(.9),
    ℓ = positive(.9)
)
kern(θ) = θ.σ^2 * with_lengthscale(SqExponentialKernel(), θ.ℓ) #Squared exponential kernel (note that hyperparameters are on the log scale)

# noise
noise_params = (;noise = fixed(exp(-2.)))   

θ = (;mean = mean_params, kernel = kernel_params, noise = noise_params)
flat_params, unflatten = value_flatten(θ)

random_inputs = filter(i -> isa(i, RandomUQInput), [x])
random_names = names(random_inputs)

gpr = gaussianprocess(
    [x],
    y,
    :y,
    kern,
    mZero,
    kernel_params,
    mean_params,
    noise_params,
    exp_design,
    LBFGS()
)

y_gp = gpr.gp(sample(x, 10)[:, 1])

function plotdata()
    plot(; xlabel="x", ylabel="y", legend=:bottomright)
    return scatter!(df[:, 1], df[:, 2]; label="training data", ms=2, markerstrokewidth=0)
end

plot_gp!(f; label) = plot!(f(sort!(sample(x, 100)[:, 1])); ribbon_scale=2, linewidth=1, label)

plotdata()
plot_gp!(gpr.gp; label="posterior f(⋅)")