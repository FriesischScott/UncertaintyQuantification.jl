using UncertaintyQuantification
using DataFrames

prior(df) = pdf(Normal(), df.x)
likelihood(df) = pdf(Normal(2, 0.5), df.x)

proposal = Normal()
x0 = (x=0.0,)
n = 10000
burnin = 100

mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

mh_samples, α = bayesianupdating(prior, likelihood, mh)

@show mean(mh_samples.x), std(mh_samples.x), α
