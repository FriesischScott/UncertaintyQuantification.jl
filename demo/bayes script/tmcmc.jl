using LinearAlgebra, DataFrames, FiniteDifferences, Dierckx, Reexport,
Distributions, Random, Plots, StatsBase

function MH(
    likelihood::Function,
    prior::Function,
    propdistsample::Function,
    propdistpdf::Function,
    startval::Array,
    mcit::Int64,
    burn::Int64
    )
    sample = zeros(mcit,size(startval, 1))
    sample[1,:] = startval

    for i in 1:mcit-1
        curr = sample[i,:]

        proposalval = curr + propdistsample() # hier curr input funktion?

        A = (likelihood(proposalval)/likelihood(curr)*prior(proposalval)/
            prior(curr)*propdistpdf(proposalval-curr)/
            propdistpdf(curr-proposalval))[1]


        if A>=rand(Uniform(0,1))
            sample[i+1,:] = proposalval
        else
            sample[i+1,:] = curr
        end
    end
    return sample[burn+1:mcit,:]
end

function resample(sset, hatw)
    snum = size(sset,1)
    rsset = zeros(snum, size(sset,2))
    hatwn = zeros(snum, 1)
    number = [1:snum...]

    for i in 1:snum
        resample = sample(number, ProbabilityWeights(vec(hatw)))
        rsset[i,:] = sset[resample, :]
        hatwn[i] = hatw[resample]
    end
    return rsset, hatwn
end

function calcbeta(
    likelihood::Function,
    betaj::Float64,
    sset::Array,
    snum::Int64
    )

    betaj1 = betaj + (2 - betaj)*0.5
    likelihoodcov = zeros(snum, 1)
    step = (2 - betaj)*0.25
    cov = 0

    for i in 1:snum
        likelihoodcov[i] = likelihood(sset[i,:])
    end

    while cov> 1.0001 || cov < 0.9999
        cov = std(likelihoodcov.^(betaj1-betaj))/mean(likelihoodcov.^(betaj1-betaj))
        if cov < 1
            betaj1 += step
            step = 0.5*step
        else
            betaj1 -= step
            step = 0.5*step
        end
        if betaj1>1
            break
        end
    end
    betaj1 = min(1, betaj1)
    return betaj1, likelihoodcov.^(betaj1-betaj)
end

function calcsigma(
    sset,
    hatw,
    gamma
    )
    sdim = size(sset, 2)
    snum = size(sset, 1)

    thet = zeros(sdim, 1)
    sig = zeros(sdim, sdim)

    for i in 1:snum
        thet = thet + sset[i,:].*hatw[i]
    end

    for i in 1:snum
        sig = sig + ((sset[i,:]-thet)*(sset[i,:]-thet)').*hatw[i]
    end
    sig = sig.*gamma^2

    return sig
end

function tmcmc(
    likelihood::Function,
    prior::Function,
    priorsample::Function,
    snum::Int64,
    gamma::Float64
    )

    j = 0
    betaj = 0.0
    sset = priorsample(snum)


    while betaj < 1
        betaj1, w = calcbeta(likelihood, betaj, sset, snum)
        hatw = w./sum(w)

        sset, w = resample(sset, hatw)
        hatw = w./sum(w)

        sig = calcsigma(sset, hatw, gamma)

        for i in 1:snum
            mvpropdist = MvNormal(zeros(size(sset, 2)), sig)
            propdistsample() = rand(mvpropdist)
            propdistpdf(x) = pdf(mvpropdist, x)
            sset[i,:] = MH(likelihood, prior, propdistsample, propdistpdf, sset[i,:], 2,1)
        end
        betaj = betaj1
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
gamma = 0.2

#rsample = tmcmc(expfunc, prior, priorsample, mcit, gamma)
#display(histogram(rsample[:, 1]))


mu = [0, 0]
sig = [1 0.8; 0.8 1400]
mvnormal = MvNormal(mu, sig)
mvnormal2 = MvNormal([3, -5], sig)
likenormal(x) = pdf(mvnormal, [x[1], x[2]])+pdf(mvnormal2, [x[1], x[2]])

mpriorsample(n) = rand(Uniform(-400,400),(n, 2))

tmsample = tmcmc(likenormal, prior, mpriorsample, 10000, gamma)

#display(histogram(msample[:, 1]))
display(scatter(tmsample[1:10000,1], tmsample[1:10000, 2]))
