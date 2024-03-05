using UncertaintyQuantification, DataFrames

l = Parameter(1.8, :l) # length
b = Parameter(0.12, :b) # width

h = RandomVariable(Normal(0.24, 0.01), :h) # height

μ, σ = distribution_parameters(10e9, 1.6e9, LogNormal)
E = RandomVariable(LogNormal(μ, σ), :E) # young's modulus

# μ, σ = distribution_parameters(5000, 400, LogNormal)
# P = RandomVariable(LogNormal(μ, σ), :P) # tip load
P = Parameter(5000, :P)
# μ, σ = distribution_parameters(600, 140, LogNormal)
# ρ = RandomVariable(LogNormal(μ, σ), :ρ) # density
ρ = Parameter(600, :ρ)

# c = GaussianCopula([1 0.8; 0.8 1])
# jd = JointDistribution([E, ρ], c)

inputs = [l, b, h, P, E, ρ]

inertia = Model(df -> df.b .* df.h .^ 3 / 12, :I)

displacement = Model(
    df ->
        (df.ρ .* 9.81 .* df.b .* df.h .* df.l .^ 4) ./ (8 .* df.E .* df.I) .+
        (df.P .* df.l .^ 3) ./ (3 .* df.E .* df.I),
    :w,
)

max_displacement = 0.01

# Compute probability of failure using standard Monte Carlo
mc = (10^6)

mc_pf, mc_std, mc_samples = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, mc
)

println(
    "Monte Carlo probability of failure $mc_pf ($(size(mc_samples, 1)) model evaluations)"
)

# Compute probability of failure using Importance Sampling
pf, β, dp, α = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, FORM()
)
is = ImportanceSampling(10^4, β, dp, α)

is_pf, is_std, is_samples = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, is
)

println(
    "Importance Sampling probability of failure: $is_pf ($(size(is_samples, 1)) model evaluations)",
)

# Compute probability of failure using Line Sampling
ls = LineSampling(200)

ls_pf, ls_std, ls_samples = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, ls
)

println(
    "Line Sampling probability of failure: $ls_pf ($(size(ls_samples, 1)) model evaluations)",
)

# Compute probability of failure using Subset Sampling
subset = UncertaintyQuantification.SubSetSimulation(2000, 0.1, 10, Uniform(-0.5, 0.5))

subset_pf, subset_std, subset_samples = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, subset
)

println(
    "Subset Simulation probability of failure: $subset_pf ($(size(subset_samples, 1)) model evaluations)",
)

# Compute probability of failure using conditional Subset Sampling
subset_inf = UncertaintyQuantification.SubSetInfinity(2000, 0.1, 10, 0.5)

subset_pf_inf, subset_std_inf, subset_samples_inf = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, subset_inf
)

println(
    "Subset infinity probability of failure: $subset_pf_inf ($(size(subset_samples_inf, 1)) model evaluations)",
)

# Compute probability of failure using adaptive Subset Sampling
subset_adap = UncertaintyQuantification.SubSetInfinityAdaptive(2000, 0.1, 10, 10, 0.6, 1.0)

subset_pf_adap, subset_std_adap, subset_samples_adap = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, subset_adap
)

println(
    "Subset infinity adaptive probability of failure: $subset_pf_adap ($(size(subset_samples_adap, 1)) model evaluations)",
)


## In SNS space
dp_vals = values(dp)
dp_vars = ["h", "E"]

dp_data_frame = DataFrame([:h => dp.h, :E => dp.E,]) 
to_standard_normal_space!(inputs, dp_data_frame)

R2 = sum(Matrix(dp_data_frame).^2)
pf_upper = 1 - cdf(Chisq(2), R2)

println("Chisquared trick pf <= $pf_upper")

in_circle(df) = df.h .^2 + df.E .^2 .<= R2
is_in_circle = in_circle(samples)

# ## In physical space
# dp_vals = values(dp)
# dp_data_frame = DataFrame([:h => dp.h, :E => dp.E]) 
# R2 = sum(Matrix(dp_data_frame).^2)

# mc_pf, mc_std, mc_samples = probability_of_failure(
#     [inertia, displacement], df -> df.h .^2 + df.E .^2 .- R2, inputs, mc
# )

# get_dist(df) =  df.h .^2 + df.P .^2 + df.E .^2 + df.ρ .^2
# R_dist = get_dist(mc_samples)

### Make plots

limit_state(df) = max_displacement .- df.w
samples = sample(inputs, SobolSampling(10^4))
evaluate!([inertia, displacement], samples)

g = limit_state(samples)
failures = g .<= 0

to_standard_normal_space!(inputs, samples)

in_circle(df) = df.h .^2 + df.E .^2 .<= R2
is_in_circle = in_circle(samples)

using Plots

SNS_centre = DataFrame([:h => 0, :E => 0])
to_physical_space!([h, E], SNS_centre)
to_physical_space!(inputs, samples)

h_samples = samples.h
E_samples = samples.E

scatter(h_samples[ .!failures], E_samples[ .!failures], color = "blue")
scatter!(h_samples[failures], E_samples[failures], color = "red")
scatter!(h_samples[is_in_circle], E_samples[is_in_circle], color = "purple")
scatter!(h_samples[h_in_c], E_samples[h_in_c], color = "green")
# scatter!([dp.h], [dp.E], color = "yellow", markershape = :star, markersize = 10)
# scatter!([SNS_centre.h], [SNS_centre.E], color = "green", markershape = :star, markersize = 10)


## In SNS space
to_standard_normal_space!(inputs, samples)
# to_physical_space!([h, E], SNS_centre)

h_samples = samples.h
E_samples = samples.E

scatter(h_samples[ .!failures], E_samples[ .!failures], color = "blue")
scatter!(h_samples[failures], E_samples[failures], color = "red")
scatter!(h_samples[is_in_circle], E_samples[is_in_circle], color = "purple")

## In SNS space
dp_vals = values(dp)
dp_vars = ["h", "E"]

dp_data_frame = DataFrame([:h => dp.h, :E => dp.E,]) 
to_standard_normal_space!(inputs, dp_data_frame)

scatter!(dp_data_frame.h, dp_data_frame.E, color = "yellow", markershape = :star, markersize = 10)
scatter!([0], [0], color = "green", markershape = :star, markersize = 10,  aspect_ratio=1)

## Copula space

h_samples_copula = cdf(Normal(), samples.h)
E_samples_copula = cdf(Normal(), samples.E)

dp_h_c = cdf(Normal(), dp_data_frame.h)
dp_E_c = cdf(Normal(), dp_data_frame.E)

scatter(h_samples_copula[ .!failures], E_samples_copula[ .!failures], color = "blue")
scatter!(h_samples_copula[failures], E_samples_copula[failures], color = "red")
scatter!(h_samples_copula[is_in_circle], E_samples_copula[is_in_circle], color = "purple")

dists_cop(x, y)  = (x .- 0.5) .^2 .+ (y .- 0.5) .^2
R2_cop = dists_cop(dp_h_c, dp_E_c)

R2_failures = dists_cop(h_samples_copula[failures], E_samples_copula[failures])
dist, index = findmin(R2_failures)
dp_copula = [h_samples_copula[failures][index], E_samples_copula[failures][index]]

C_samples = rand(10^5, 2)

C_dists = dists_cop(C_samples[:,1], C_samples[:, 2])
C_in = C_dists .<= dist
h_in_c = dists_cop(h_samples_copula, E_samples_copula) .<= dist

scatter!(C_samples[C_in, 1], C_samples[C_in, 2], color = "green")
scatter!([dp_copula[1]], [dp_copula[2]], color = "yellow", markershape = :star, markersize = 10,  aspect_ratio=1)