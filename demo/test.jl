using UncertaintyQuantification

function expfunc(x)
    if x[1]<0
        rv = 0
    else
        rv = exp(-x[1])
    end
    return rv
end
prior(x) = 1
propdist = Normal(0, 1)
propdistsample() = [rand(propdist)]
propdistpdf(x) = pdf(propdist, x[1])



mcit = MonteCarlo(50500)

rsample = MH(expfunc, prior, propdistpdf, propdistsample, [3], mcit, 500)

display(histogram(rsample[:, 1]))