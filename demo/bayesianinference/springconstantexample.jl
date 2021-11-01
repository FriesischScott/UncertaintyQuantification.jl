using UncertaintyQuantification

Random.seed!(4183);

N = 10
normal = Normal(0, 0.005)

F = [5, 5, 10, 10, 15, 15, 20, 20, 25, 25]./1.02

noisesamples = [F[i]/387 + rand(normal) for i in 1:N]

observations = hcat(F,noisesamples)
Random.seed!()

function likelihood(x)
    rv = 0
    variation = var(observations[:,1]./observations[:,2])
    for i in 1:N
        rv += (observations[i,1]/observations[i,2] - x[1])^2
    end
    return exp(-rv/(2*variation))
 end

prior(x) = 1
priorsample(n) = rand(Uniform(1,1000),(n))

propdist = Normal(0, 1)
propdistsample() = [rand(propdist)]
propdistpdf(x) = pdf(propdist, x[1])

function term(sset, j)
    rv = true
    if j >= 2 rv = false end
    return rv
end


smcsample = smc(likelihood, prior, priorsample, propdistpdf, propdistsample, 5000, term)