###
#   Tests accuracy of MCMC sampler vs analytical solution provided by conjugate priors: https://en.wikipedia.org/wiki/Conjugate_prior
###

using UncertaintyQuantification
using Plots

mean_true = 5
std_fixed = 2     # Fixed

N_data = 5

Data = rand(Normal(mean_true, std_fixed), N_data)

prior_mean = 2
prior_std = 10

# Exact conjugate posterior
posterior_std = 1/((1/prior_std^2) + (N_data/std_fixed^2))
posterior_mean = posterior_std * (prior_mean / prior_std^2 + sum(Data) /std_fixed^2)
posterior_exact(x) = pdf(Normal(posterior_mean, posterior_std), x)


prior(df) = pdf(Normal(prior_mean, prior_std), df.x)
likelihood(df) = prod(pdf.(Normal.(df.x, std_fixed), Data))


# Proposal
proposal = Normal(0, 1.5)

x0 = (x=2.0,)
n = 20000
burnin = 200

mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)
mh_samples, α = bayesianupdating(prior, likelihood, mh)

xs = range(-1, 10, length = 100)
histogram(mh_samples.x, normalize=:pdf, bin = 100, label = "MCMC")
plot!(xs, posterior_exact.(xs), linewidth = 2, label = "Congugate prior", xlabel = "normal mean")
vline!([mean_true], label = "mean true", linewidth = 3, color = "black")

println("Exact mean: $posterior_mean | exact std: $posterior_std")
println("TMCMC mean: $(mean(mh_samples.x)) | TMCMC std: $(std(mh_samples.x))")