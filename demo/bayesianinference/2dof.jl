using UncertaintyQuantification
using Plots

λ1 = @. Model(df -> ((df.θ1 + 2 * df.θ2) + sqrt(df.θ1^2 + 4 * df.θ2^2)) / 2, :λ1)
λ2 = @. Model(df -> ((df.θ1 + 2 * df.θ2) - sqrt(df.θ1^2 + 4 * df.θ2^2)) / 2, :λ2)

function likelihood(df, data::AbstractMatrix)
    σ = [1 0.1]
    λ = [df.λ1 df.λ2]
    return exp(-0.5 * sum(((data .- λ) ./ σ) .^ 2))
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
    return pdf(Uniform(0, 4), df.θ1) .* pdf(Uniform(0, 4), df.θ2)
end

lh(df) = likelihood(df, data)

proposal = Normal(0, 0.1)
x0 = (θ1=2.84, θ2=2.33)
n = 10000
burnin = 450

mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

mh_samples, α = bayesianupdating(prior, Like, [λ1, λ2], mh)

println("Acceptance rate: $α")
println("Identified (θ1, θ2): ($(mean(mh_samples.θ1)), $(mean(mh_samples.θ2)))")

scatter(
    mh_samples.θ1,
    mh_samples.θ2;
    xlabel="θ1",
    ylabel="θ2",
    color="red",
    label="Posterior Input",
)
