using UncertaintyQuantification
using Plots
using Random

Random.seed!(2873)

prior_sample = RandomVariable.(Uniform(-5, 5), [:x, :y])

prior = df -> logpdf.(Uniform(-5, 5), df.x) .+ logpdf(Uniform(-5, 5), df.y)

himmelblau = Model(
    df -> (df.x .^ 2 .+ df.y .- 11) .^ 2 .+ (df.x .+ df.y .^ 2 .- 7) .^ 2, :f
)

likelihood = df -> -df.f

n = 2000

burnin = 20

tmcmc = UncertaintyQuantification.TMCMC(prior_sample, n, burnin, 0.2)

samples, S = bayesianupdating(prior, likelihood, [himmelblau], tmcmc)

scatter(samples.x, samples.y; aspect_ratio=:equal, lims=[-5, 5])

@show S
