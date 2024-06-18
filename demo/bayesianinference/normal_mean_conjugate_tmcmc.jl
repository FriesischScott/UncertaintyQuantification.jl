###
#   Tests accuracy of TMCMC sampler vs analytical solution provided by conjugate priors: https://en.wikipedia.org/wiki/Conjugate_prior
###

using UncertaintyQuantification
using Plots

mean_true = 4
std_fixed = 5     # Fixed

N_data = 25

Data = rand(Normal(mean_true, std_fixed), N_data)

prior_mean = 2
prior_std = 10

# Exact conjugate posterior
posterior_std = 1/((1/prior_std^2) + (N_data/std_fixed^2))
posterior_mean = posterior_std * (prior_mean / prior_std^2 + sum(Data) /std_fixed^2)
posterior_exact(x) = pdf(Normal(posterior_mean, posterior_std), x)


# Proposal
proposal = Normal(0, 1)

# TMCMC 
logprior(df) = logpdf(Normal(prior_mean, prior_std), df.x)
prior_sample_ = RandomVariable(Normal(prior_mean, prior_std), :x)
loglikelihood(df) = [sum(logpdf.(Normal.(df_i.x, std_fixed), Data)) for df_i in eachrow(df)]


n = 2000
burnin = 200
thin = 5

tmcmc = UncertaintyQuantification.TMCMC([prior_sample_], n, burnin, thin, 0.01)
samples, log_evidence = bayesianupdating(logprior, loglikelihood, UQModel[], tmcmc)

xs = range(-1, 10, length = 100)
histogram(samples.x, normalize=:pdf, bin = 100, label = "MCMC")
plot!(xs, posterior_exact.(xs), linewidth = 2, label = "Congugate prior", xlabel = "normal mean")
vline!([mean_true], label = "mean true", linewidth = 3, color = "black")

println("Exact mean: $posterior_mean | exact std: $posterior_std")
println("TMCMC mean: $(mean(samples.x)) | TMCMC std: $(std(samples.x))")