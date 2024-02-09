using UncertaintyQuantification
using Plots
using LaTeXStrings
using DataFrames
using StatsPlots

function analyticalexample_2D_noise(θ::AbstractArray)
    λ1_model = ((θ[1] + 2 * θ[2]) + sqrt(θ[1]^2 + 4 * θ[2]^2)) / 2
    λ2_model = ((θ[1] + 2 * θ[2]) - sqrt(θ[1]^2 + 4 * θ[2]^2)) / 2

    λ1_noise = λ1_model + rand(Normal())
    λ2_noise = λ2_model + rand(Normal(0, 0.5))
    return [λ1_noise λ2_noise]
end

function analyticalexample_2D(θ::AbstractArray)
    λ1_model = ((θ[1] + 2 * θ[2]) + sqrt(θ[1]^2 + 4 * θ[2]^2)) / 2
    λ2_model = ((θ[1] + 2 * θ[2]) - sqrt(θ[1]^2 + 4 * θ[2]^2)) / 2
    return [λ1_model λ2_model]
end

function likelihood(θ::AbstractArray, model::Function, Yexp::AbstractMatrix)
    σ1 = 1
    σ2 = 0.5
    Ysim = model(θ)
    temp_exp =
        -0.5 *
        sum((((Yexp[:, 1] .- Ysim[1]) / σ1) .^ 2 + ((Yexp[:, 2] .- Ysim[2]) / σ2) .^ 2))
    return exp(temp_exp)
end

### Obervations
#N_obs = 15
#λ_obs = Matrix{Float64}(undef, N_obs, 2)
λ_obs = [
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
#for i_obs in 1:N_obs
#    λ_obs[i_obs, :] = analyticalexample_2D_noise([0.5 1.5])
#end

### Prior 
function prior(x)
    return pdf(Uniform(0, 4), x[1]) * pdf(Uniform(0, 4), x[2])
end

### True value
λ_true = analyticalexample_2D([0.5 1.5])

### Bayesian Model Updating 
Like(x) = likelihood(x, analyticalexample_2D, λ_obs)

proposal = Normal()
x0 = DataFrame(:x => 2.84, :y => 2.33)
n = 1000
burnin = 100

mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

mh_samples = bayesianupdating(prior, Like, mh)

@show mean(mh_samples.x), std(mh_samples.y)

#corrplot(Matrix(mh_samples))
### Figures
scatter(
    λ_obs[:, 1],
    λ_obs[:, 2];
    xlabel="Eigenvalue " * L"λ_1^{noisy}",
    ylabel="Eigenvalue " * L"λ_2^{noisy}",
    color="red",
    label="Noisy eigenvalue",
)

scatter!([λ_true[1, 1]], [λ_true[1, 2]]; color="black", label="True eigenvalue", marker=:+)

scatter(
    mh_samples.x,
    mh_samples.y;
    xlabel=L"θ_1",
    ylabel=L"θ_2",
    color="red",
    label="Posterior Input",
)