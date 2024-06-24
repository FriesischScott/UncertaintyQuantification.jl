using UncertaintyQuantification
using Plots
using Random

Random.seed!(2873)

prior = RandomVariable.(Uniform(-5, 5), [:x, :y])

himmelblau = Model(
    df -> (df.x .^ 2 .+ df.y .- 11) .^ 2 .+ (df.x .+ df.y .^ 2 .- 7) .^ 2, :f
)

likelihood = df -> -df.f

n = 2000

burnin = 20

tmcmc = TransitionalMarkovChainMonteCarlo(prior, n, burnin, 0.2)

samples, S = bayesianupdating(likelihood, [himmelblau], tmcmc)

@show S

scatter(samples.x, samples.y; aspect_ratio=:equal, lims=[-5, 5])
