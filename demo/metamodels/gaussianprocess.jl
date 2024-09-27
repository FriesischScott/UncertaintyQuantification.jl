using UncertaintyQuantification, Plots
using AbstractGPs # not sure if to reexport
using ParameterHandling # not sure if to reexport
using Optim
using Random

Random.seed!(20140430)
n = 10
input_symbol = :x
output_symbol = :y

# Input_symbol for passing a UQModel to GaussianProcess
x = RandomVariable(Uniform(0, 2π), input_symbol)
y = Model(
    df ->
        (sin.(df.x) + 0.05*randn(length(df.x))),
        output_symbol,
)
exp_design = ExperimentalDesign(MonteCarlo(n))

# Input_symbol for passing a DataFrame to GaussianProcess
df = sample(x, n)
evaluate!(y, df)

# Define how the GP is build
function build_gp(θ::NamedTuple)
    k1 = (θ.SE.σ)^2 * with_lengthscale(SqExponentialKernel(), θ.SE.ℓ)
    k2 = (θ.RQ.σ)^2 * with_lengthscale(RationalQuadraticKernel(; α=θ.RQ.α), θ.RQ.ℓ)
    return GP(k1+k2)
end

# gp parameters
θ = (;
    SE = (;
        σ = positive(1.),
        ℓ = positive(1.)
    ),
    RQ = (;
        σ = positive(1.),
        ℓ = positive(1.),
        α = positive(exp(-1.0))
    )
)

# noise
noise = fixed(exp(-2.)) 

# Fit GP to data from experimental design
gp_from_model = GaussianProcess(
    x,
    y,
    :y,
    exp_design,
    build_gp,
    θ,
    noise
)

# Optimize hyperparameters (There should be a method that allows to do this on a already fitted gp instance)
opt_gp_from_model = GaussianProcess(
    x,
    y,
    :y,
    exp_design,
    build_gp,
    θ,
    noise,
    false,
    false,
    LBFGS()
)

function plotdata()
    plot(; xlabel="x", ylabel="y", legend=:bottomright)
    return scatter!(df[:, 1], df[:, 2]; label="training data", ms=2, markerstrokewidth=0)
end

plot_gp!(f; label) = plot!(f(sort!(sample(x, 100)[:, 1])); ribbon_scale=2, linewidth=1, label)

plotdata()
plot_gp!(gp_from_model.gp; label="posterior f(⋅)")
plot_gp!(opt_gp_from_model.gp; label="posterior f(⋅) optimized")