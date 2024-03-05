using UncertaintyQuantification, Plots

# u1 = RandomVariable(Normal(0, 1), :u1)
# u2 = RandomVariable(Normal(0, 1), :u2)

u1 = RandomVariable(Uniform(-10, 10), :u1)
u2 = RandomVariable(Uniform(-10, 10), :u2)

inputs = [u1, u2]

function performance(df)
    u1 = df.u1
    u2 = df.u2

    p1 = 0.1 .* (u1 .- u2) .^2 .- (u1 .+ u2) ./sqrt(2) .+3
    p2 = 0.1 .* (u1 .- u2) .^2 .+ (u1 .+ u2) ./sqrt(2) .+3
    p3 = u1 .- u2 .+ 7*sqrt(2)
    p4 = u2 .- u1 .+ 7*sqrt(2)

    return minimum([p1 p2 p3 p4], dims = 2)
end

# displacement = Model(
#     df ->
#         (df.ρ .* 9.81 .* df.b .* df.h .* df.l .^ 4) ./ (8 .* df.E .* df.I) .+
#         (df.P .* df.l .^ 3) ./ (3 .* df.E .* df.I),
#     :w,
# )

# max_displacement = 0.01

# Compute probability of failure using standard Monte Carlo
mc = SobolSampling(10^4)

mc_pf, mc_std, mc_samples = probability_of_failure( performance, inputs, mc)

pf, β, dp, α= probability_of_failure( performance, inputs, FORM())

failures = performance(mc_samples) .<= 0
failures = reduce(vcat, failures)

scatter(mc_samples.u1[ .!failures], mc_samples.u2[ .!failures], color= "blue")
scatter!(mc_samples.u1[failures], mc_samples.u2[failures], color= "red")

samples = mc_samples

## Functionise it

function FORM_norm(samples, failures, p)

    norms = [norm(samples[i,:], p) for i = 1:size(samples, 1)]

    norms_failures = norms[failures]
    failure_samples = samples[failures, :]

    distance, dp_index = findmin(norms_failures)

    inside_index = norms .<= distance
    
    return inside_index, failure_samples[dp_index,:]
end

samples_SNS = Matrix(samples[!, [:u1, :u2]])


plt = scatter(samples_SNS[ .!failures, 1], samples_SNS[ .!failures, 2], color = "blue")
scatter!(plt, samples_SNS[failures, 1], samples_SNS[failures, 2], color = "red")

ps = [0.01, 0.1, 0.5, 1, 1.3, 1.5, 1.7, 2, Inf]
for p in ps 
    inside_index, dp = FORM_norm(samples_SNS, failures, p)

    scatter!(plt, samples_SNS[inside_index, 1], samples_SNS[inside_index, 2], label = "p = $p")
    scatter!(plt, [dp[1]], [dp[2]], label = "dp = $p", markershape = :star, markersize=14)

end

plt



function l2_to_my_norm(samples, p1, p2)

    Ns = size(samples, 1)
    norm_p1 = [norm(samples[i,:], p1)  for i in 1:Ns] 
    norm_p1 = reduce(vcat, norm_p1)

    samples_normlised = samples ./ norm_p1

    norm_p2 = [norm(samples[i,:], p2)  for i in 1:Ns] 

    return samples_normlised .* norm_p2
end

samples_SNS_transformed = lp_to_lp(samples_SNS, 2, 0.2)
# samples_SNS_transformed_failures = samples_SNS_transformed[failures]

plt = scatter(samples_SNS_transformed[ .!failures, 1], samples_SNS_transformed[ .!failures, 2], color = "blue")
scatter!(plt, samples_SNS_transformed[failures, 1], samples_SNS_transformed[failures, 2], color = "red", aspect_ratio=1)


