using UncertaintyQuantification, Plots

# Model parameters
L1 = Parameter(3, :L1)
L2 = Parameter(2, :L2)
Zp = Parameter(3.66e-4, :Zp)

# Random variables
β = 5 * sqrt(6) / pi;
μ = 50 - 0.5772 * β;
X1 = RandomVariable(Gumbel(μ, β), :X1)

# Yield stress
lb = 19.9e4 # lower bound
μ_phys = 28.8e4 # physical mean μ
σ_phys = 2.64e4 # physical standard deviation σ

μ = log((μ_phys)^2 / (sqrt(σ_phys^2 + (μ_phys)^2)))
σ = sqrt(log(σ_phys^2 / (μ_phys)^2 + 1))
X2 = RandomVariable(Truncated(LogNormal(μ, σ), lb, Inf), :X2)

inputs = [L1, L2, Zp, X1, X2]

performance = Model(df -> (-df.L1 .* df.L2) ./ (df.L1 + df.L2) .* (1 ./ df.Zp) .* df.X1 .+ df.X2, :p)

# Compute probability of failure using standard Monte Carlo
mc = MonteCarlo(10^5)

mc_pf, mc_samples = probability_of_failure(
    [performance],
    df -> df.p,
    inputs,
    mc,
)

sus = SubSetSimulation(2000, 0.2, 10, Uniform(-0.5, 0.5))

sus_pf, sus_samples = probability_of_failure(
    [performance],
    df -> df.p,
    inputs,
    sus,
)

p = plot(aspect_ratio=:equal)

for level ∈ unique(sus_samples.level)
    data = filter(:level => l -> l == level, sus_samples)
    to_standard_normal_space!([X1, X2], data)
    scatter!(p, data.X1, data.X2, label="l$level")
end

display(p)