using UncertaintyQuantification

l = Parameter(1.8, :l) # length
b = Parameter(0.12, :b) # width

h = RandomVariable(Normal(0.24, 0.01), :h) # height

μ = log(10e9^2 / sqrt(1.6e9^2 + 10e9^2))
σ = sqrt(log(1.6e9^2 / 10e9^2 + 1))
E = RandomVariable(LogNormal(μ, σ), :E) # young's modulus

μ = log(5000^2 / sqrt(400^2 + 5000^2))
σ = sqrt(log(400^2 / 5000^2 + 1))
P = RandomVariable(LogNormal(μ, σ), :P) # tip load

μ = log(600^2 / sqrt(140^2 + 600^2))
σ = sqrt(log(140^2 / 600^2 + 1))
ρ = RandomVariable(LogNormal(μ, σ), :ρ) # density

c = GaussianCopula([1 0.8; 0.8 1])
jd = JointDistribution([E ρ], c)

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

subset_pf, subset_cov, subset_samples= probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, subset
)

println(
    "Subset Simulation probability of failure: $subset_pf ($(size(subset_samples, 1)) model evaluations)",
)
