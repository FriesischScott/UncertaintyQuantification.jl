using UncertaintyQuantification
using AbstractGPs
using Random

# Setup Himmelblau example
x = RandomVariable.(Uniform(-5, 5), [:x1, :x2])
himmelblau = Model(
    df -> (df.x1 .^ 2 .+ df.x2 .- 11) .^ 2 .+ (df.x1 .+ df.x2 .^ 2 .- 7) .^ 2, :y
)
design = FullFactorial([8, 8])
training_data = sample(x, design)
evaluate!(himmelblau, training_data)

# This will be used for proper initial guesses for the parameters of the GP
mean_data = mean(training_data[!, :y])
std_data = std(training_data[!, :y])

# Setup the GP
# Note: If we do not initialize the parameters here properly the optimization will fail. Standardization should help with that.
σ² = 1e-5
kernel = SqExponentialKernel() ∘ ARDTransform([1.0, 1.0])
gp = with_gaussian_noise(GP(0.0, kernel), σ²)

optimizer = MLE()
# TODO: StandardizeInput breaks currently due to -Inf and Inf values from to_standard_normal_space!()
# TODO: Optimization is extremely unstable
# TODO: Not all kernels have a extract_parameters and apply_parameters function 
gpr = GaussianProcess(
    gp, x, himmelblau, :y, design, StandardizeOutput(), MLE()
)

test_data = sample(x, 1000)
evaluate!(gpr, test_data)

p_data = test_data[:, [:x1, :x2]]
evaluate!(himmelblau, p_data)

mse = mean((p_data.y .- test_data.y) .^ 2)
println("MSE is:  $mse")

using Plots
using DataFrames
a = range(-5, 5; length=1000) 
b = range(5, -5; length=1000)
himmelblau_f(x1, x2) = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
function gpr_f(x, y)
    df = DataFrame(x1 = x, x2 = y)
    evaluate!(gpr, df)
    return only(df[:, :y])
end

s1 = surface(a, b, himmelblau_f; plot_title="Himmelblau's function")
s2 = surface(a, b, gpr_f; plot_title="Gaussian process regression")
plot(s1, s2, layout = (1, 2), legend = false)