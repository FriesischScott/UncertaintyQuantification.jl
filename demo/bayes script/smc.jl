using LinearAlgebra, DataFrames, FiniteDifferences, Dierckx, Reexport,
Distributions, Random, Plots, StatsBase

function resample(sset, hatw)
    snum = size(sset,1)
    rsset = zeros(snum, size(sset,2))
    number = [1:snum...]

    for i in 1:snum
        resample = sample(number, ProbabilityWeights(vec(hatw)))
        rsset[i,:] = sset[resample, :]
    end

    return rsset
end

function smc(
    likelihood::Function,
    prior::Function,
    priorsample::Function,
    propdistpdf::Function,
    propdistsample::Function,
    snum::Int64,
    term::Float64
    )

    j = 0
    cov = 0
    sset = priorsample(snum)
    w = zeros(snum,1)

    for i in 1:snum
        w[i] = likelihood(sset[i,:])*prior(sset[i,:])
    end

    hatw = w./sum(w)

    while cov < term
        j += 1
        Neff = 1/sum(hatw.^2)
        if Neff < snum/2
            sset = resample(sset, hatw)
            hatw = ones(snum,1)./snum
        end

        rsset = zeros(snum, size(sset,2))
        w = zeros(snum, 1)
        likelihoodcov = zeros(snum, size(sset,2))

        for i in 1:snum
            rsset[i,:] = sset[i,:] + propdistsample()
            likelihoodcov[i] = likelihood(rsset[i,:])/likelihood(rsset[i,:])
            w[i] = hatw[i].*likelihood(rsset[i,:])/likelihood(sset[i,:])*
                prior(rsset[i,:])/prior(sset[i,:])*
                propdistpdf(rsset[i,:]-sset[i,:])/propdistpdf(sset[i,:]-rsset[i,:])
        end

        sset = rsset
        hatw = w./sum(w)

        cov = std(likelihoodcov)/mean(likelihoodcov)

    end
    return sset
end


function expfunc(x)
    if x[1]<0
        rv = 0
    else
        rv = exp(-x[1])
    end
    return rv
end
prior(x) = 1

priorsample(n) = rand(Uniform(0,10),n)

propdist = Normal(0, 0.1)
propdistsample() = [rand(propdist)]
propdistpdf(x) = pdf(propdist, x[1])


mcit = 5000
term = 0.2

#rsample = smc(expfunc, prior, priorsample, propdistpdf, propdistsample, mcit, term)

#display(histogram(rsample[:, 1]))


mu = [0, 0]
sig = [1 0; 0 20]
mvnormal = MvNormal(mu, sig)
mvnormal2 = MvNormal([11, -7], sig)
likenormal(x) = pdf(mvnormal, [x[1], x[2]]) +pdf(mvnormal2, [x[1], x[2]])

mpriorsample(n) = rand(Uniform(-30,30),(n, 2))


propsig = [0.5 0; 0 0.5]
mpropdist = MvNormal(mu, propsig)
mpropdistsample() = rand(mpropdist)
mpropdistpdf(x) = pdf(mpropdist, x)
smsample = smc(likenormal, prior, mpriorsample, mpropdistpdf, mpropdistsample, 50000, term)

#display(histogram(msample[:, 1]))
display(scatter(smsample[1:10000,1], smsample[1:10000, 2]))
#eqt = calcequaltails(smsample, 0.8)
