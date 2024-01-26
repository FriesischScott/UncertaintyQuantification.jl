using UncertaintyQuantification

prior(x) = pdf(Normal(), x)
likelihood(x) = pdf(Normal(2, 0.5), x)

prop_pdf(x) = pdf(Normal(), x)
prop_sample() = rand(Normal())
x0 = DataFrame(:x => 0.0, :y => 0.0)
n = 1000
burnin = 100

mh = MetropolisHastings(
    prop_pdf, prop_sample, x0, n, burnin
)

mh_samples = bayesianupdating(prior, likelihood, mh)

@show mean(mh_samples.x), std(mh_samples.x)
