using UncertaintyQuantification
using Plots
using LaTeXStrings
using StatsPlots

λ1 = @. Model(df -> ((df.θ1 + 2 * df.θ2) + sqrt(df.θ1^2 + 4 * df.θ2^2)) / 2, :λ1)
λ2 = @. Model(df -> ((df.θ1 + 2 * df.θ2) - sqrt(df.θ1^2 + 4 * df.θ2^2)) / 2, :λ2)

function likelihood(df, data::AbstractMatrix)
    σ1 = 1
    σ2 = 0.1
    temp_exp =
        -0.5 * sum((((data[:, 1] .- df.λ1) / σ1) .^ 2 + ((data[:, 2] .- df.λ2) / σ2) .^ 2))
    return exp.(temp_exp)
end

### Obervations
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

### Prior
function prior(df)
    return pdf(Uniform(0, 4), df.θ1) .* pdf(Uniform(0, 4), df.θ2)
end

### Bayesian Model Updating
Like(df) = likelihood(df, data)

proposal = Normal(0, 0.1)
x0 = (θ1=2.84, θ2=2.33)
n = 10000
burnin = 450

mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

mh_samples, α = bayesianupdating(prior, Like, [λ1, λ2], mh)

@show mean(mh_samples.θ1), std(mh_samples.θ2), α

scatter(
    mh_samples.θ1,
    mh_samples.θ2;
    xlabel=L"θ_1",
    ylabel=L"θ_2",
    color="red",
    label="Posterior Input",
)
