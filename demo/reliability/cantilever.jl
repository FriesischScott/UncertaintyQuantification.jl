using UncertaintyQuantification

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

max_displacement = 0.01

# Compute probability of failure using standard Monte Carlo
mc = MonteCarlo(10^6)

mc_pf, mc_cov, mc_samples = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, mc
)

println(
    "Monte Carlo probability of failure $mc_pf ($(size(mc_samples, 1)) model evaluations)"
)

# Compute probability of failure using Line Sampling
ls = LineSampling(200)

ls_pf, ls_cov, ls_samples = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, ls
)

println(
    "Line Sampling probability of failure: $ls_pf ($(size(ls_samples, 1)) model evaluations)",
)

subset = UncertaintyQuantification.SubSetSimulation(2000, 0.1, 10, Uniform(-0.5, 0.5))

subset_pf, subset_cov, subset_samples = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, subset
)

println(
    "Subset Simulation probability of failure: $subset_pf ($(size(subset_samples, 1)) model evaluations)",
)
