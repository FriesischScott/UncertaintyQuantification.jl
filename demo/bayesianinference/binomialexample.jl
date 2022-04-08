using UncertaintyQuantification

MHexlikelihood(x) = x[1]^85 * (1-x[1])^(100-85)

prior(x) = 1

propdist = Normal(0, 0.2)
propdistsample() = [rand(propdist)]
propdistpdf(x) = pdf(propdist, x[1])

mcit=10000

rsample = mh(MHexlikelihood, prior, propdistpdf, propdistsample, [1], mcit, 100)