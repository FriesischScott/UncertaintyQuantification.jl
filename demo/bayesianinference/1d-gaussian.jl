using UncertaintyQuantification
using DataFrames

prior(x) = pdf(Normal(), x)
likelihood(x) = pdf(Normal(2, 0.5), x)

proposal = Normal()
x0 = (x=0.0,)
n = 10000
burnin = 100

mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

mh_samples = bayesianupdating(prior, likelihood, mh)

@show mean(mh_samples.x), std(mh_samples.x)
