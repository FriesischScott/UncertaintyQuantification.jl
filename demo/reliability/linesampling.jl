using UncertaintyQuantification

# Example from:
# Dang, C., Valdebenito, M. A., Song, J., Wei, P., & Beer, M. (2023).
# Estimation of small failure probabilities by partially Bayesian active learning
# line sampling: Theory and algorithm. Computer Methods in Applied Mechanics
# and Engineering, 412, 116068. https://doi.org/10.1016/j.cma.2023.116068

#   Input:
#       𝑿 = [𝑋₁, 𝑋₂] with 𝑋ᵢ ∼ 𝒩(0,1)
#   Limit-state function:
#       𝑔(𝑿) = 5 - 𝑋₂ + 0.01 𝑋₁³ + sin(𝑋₁)

# Input
x = RandomVariable.(Normal(), [:x1, :x2])

# Model
y = Model(df -> df.x2 .+ 0.01 * df.x1.^3 .+ sin.(df.x1), :y)

# Limit-state function
g(df) = 5 .- df.y

# Simulation
numberlines = 200

# Perform Line Sampling
LS = LineSampling(numberlines, collect(0.5:0.5:10))
pf_ls, σ_pf_ls, samples_ls = probability_of_failure([y], g, x, LS)
cov_ls = σ_pf_ls/pf_ls

# Perform Advanced Line Sampling
ALS = AdvancedLineSampling(numberlines, collect(0.5:0.5:10))
pf_als, σ_pf_als, samples_als = probability_of_failure([y], g, x, ALS)
cov_als = σ_pf_als/pf_als

println("Line Sampling:")
println("Pf = $pf_ls")
println("CoV = $cov_ls")
println("Model Calls = $(length(samples_ls.x1))\n")

println("Advanced Line Sampling:")
println("Pf = $pf_als")
println("CoV = $cov_als")
println("Model Calls = $(length(samples_als.x1))")
