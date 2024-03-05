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
mc = MonteCarlo(10^6)

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

### Find design points with different norms
using Plots, LinearAlgebra

metric(x) = norm(x, 2s)

pf, β, dp_2, α = probability_of_failure(
    [inertia, displacement], df -> max_displacement .- df.w, inputs, FORM(), metric
)



# limit_state(df) = max_displacement .- df.w
# samples = sample(inputs, SobolSampling(10^4))
# evaluate!([inertia, displacement], samples)

# g = limit_state(samples)
# failures = g .<= 0

# to_standard_normal_space!(inputs, samples)

# h_samples_SNS = samples.h
# E_samples_SNS = samples.E

# h_samples_SNS_f = h_samples_SNS[failures]
# E_samples_SNS_f = E_samples_SNS[failures]


# # Find points in failure domain with different norms

# failure_samples = [h_samples_SNS_f E_samples_SNS_f]
# num_failures = size(failure_samples,1)

# # ps = [1, 1.5, 2, 5, 10, 100, Inf]
# ps = [0.1, 0.5, 1, 1.3, 1.5, 1.7, 2]

# distances = zeros(length(ps))
# indecies = zeros(length(ps))

# for (j, p) in enumerate(ps)
#     norms = [norm(failure_samples[i,:], p) for i = 1:num_failures]
#     distances[j], indecies[j] = findmin(norms)
# end

# indecies = Integer.(indecies)

# scatter(h_samples_SNS[ .!failures], E_samples_SNS[ .!failures], color = "blue")
# scatter!(h_samples_SNS_f, E_samples_SNS_f, color = "red")

# # [scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, label="p = $(ps[i])") for i = 1:length(ps)]
# # [scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, makersize = 14) for i = 1:length(ps)]

# i = 1
# scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, markersize=14, label="p = $(ps[i])")
# i = 2
# scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, markersize=14, label="p = $(ps[i])")
# i = 3
# scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, markersize=14, label="p = $(ps[i])")
# i = 4
# scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, markersize=14, label="p = $(ps[i])")
# i = 5
# scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, markersize=14, label="p = $(ps[i])")
# i = 6
# scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, markersize=14, label="p = $(ps[i])")
# i = 7
# scatter!([h_samples_SNS_f[indecies[i]]], [E_samples_SNS_f[indecies[i]]], markershape = :star, markersize=14, label="p = $(ps[i])")


## Functionise it

function FORM_norm(samples, failures, p)

    norms = [norm(samples[i,:], p) for i = 1:size(samples, 1)]

    norms_failures = norms[failures]
    failure_samples = samples[failures, :]

    distance, dp_index = findmin(norms_failures)

    inside_index = norms .<= distance
    
    return inside_index, failure_samples[dp_index,:]
end

limit_state(df) = max_displacement .- df.w
samples = sample(inputs, SobolSampling(10^4))
evaluate!([inertia, displacement], samples)

g = limit_state(samples)
failures = g .<= 0

to_standard_normal_space!(inputs, samples)

samples_SNS = Matrix(samples[!, [:h, :E]])


plt = scatter(samples_SNS[ .!failures, 1], samples_SNS[ .!failures, 2], color = "blue")
scatter!(plt, samples_SNS[failures, 1], samples_SNS[failures, 2], color = "red", aspect_ratio=1)

ps = reverse([1, 1.2, 1.5, 1.7, 2])
for p in ps 
    inside_index, dp = FORM_norm(samples_SNS, failures, p)

    scatter!(plt, samples_SNS[inside_index, 1], samples_SNS[inside_index, 2], label = "p = $p")
    scatter!(plt, [dp[1]], [dp[2]], label = "dp = $p", markershape = :star, markersize=14, aspect_ratio=1)

end

plt


my_norm(X) = 0.5 * x[1]^2 + 0.9*abs(x[2])

function lp_to_lp(samples, p1, p2)

    Ns = size(samples, 1)
    norm_p1 = [norm(samples[i,:], p1)  for i in 1:Ns] 
    norm_p1 = reduce(vcat, norm_p1)

    samples_normlised = samples ./ norm_p1

    norm_p2 = [norm(samples[i,:], p2)  for i in 1:Ns] 

    return samples_normlised .* norm_p2
end

function l2_to_my_norm(samples, p1, p2)

    Ns = size(samples, 1)
    norm_p1 = [norm(samples[i,:], p1)  for i in 1:Ns] 
    norm_p1 = reduce(vcat, norm_p1)

    samples_normlised = samples ./ norm_p1

    norm_p2 = [norm(samples[i,:], p2)  for i in 1:Ns] 

    return samples_normlised .* norm_p2
end

samples_SNS_transformed = lp_to_lp(samples_SNS, 2, 1)
# samples_SNS_transformed_failures = samples_SNS_transformed[failures]

plt = scatter(samples_SNS_transformed[ .!failures, 1], samples_SNS_transformed[ .!failures, 2], color = "blue")
scatter!(plt, samples_SNS_transformed[failures, 1], samples_SNS_transformed[failures, 2], color = "red")


