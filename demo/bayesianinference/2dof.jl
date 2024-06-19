using UncertaintyQuantification
using Plots

λ1 = @. Model(df -> ((df.θ1 + 2 * df.θ2) + sqrt(df.θ1^2 + 4 * df.θ2^2)) / 2, :λ1)
λ2 = @. Model(df -> ((df.θ1 + 2 * df.θ2) - sqrt(df.θ1^2 + 4 * df.θ2^2)) / 2, :λ2)

function likelihood(df, data::AbstractMatrix)
    σ = [1 0.1]
    λ = [df.λ1 df.λ2]

    return log.(exp.([-0.5 * sum(((data .- λ[n, :]') ./ σ) .^ 2) for n in axes(λ, 1)]))
end

# Obervations
data = [
    1.51 0.33
    4.01 0.3
    3.16 0.27
    3.21 0.18
    2.19 0.33
    1.71 0.23
    2.73 0.21
    5.51 0.2
    1.95 0.11
    4.48 0.2
    1.43 0.16
    2.91 0.26
    3.81 0.23
    3.58 0.25
    2.62 0.25
]

# Prior
function prior(df)
    return logpdf(Uniform(0, 4), df.θ1) .+ logpdf(Uniform(0, 4), df.θ2)
end

prior_sample = RandomVariable.(Uniform(0, 4), [:θ1, :θ2])

lh(df) = likelihood(df, data)

n = 1000
burnin = 0

tmcmc = UncertaintyQuantification.TMCMC(prior_sample, n, burnin, 0.2)

tmcmc_samples, S = bayesianupdating(prior, lh, [λ1, λ2], tmcmc)

println("Log evidence: $S")
println("Identified (θ1, θ2): ($(mean(tmcmc_samples.θ1)), $(mean(tmcmc_samples.θ2)))")

scatter(
    tmcmc_samples.θ1,
    tmcmc_samples.θ2;
    xlabel="θ1",
    ylabel="θ2",
    color="red",
    label="",
    xlim=[0, 4],
    ylim=[0, 2],
)
