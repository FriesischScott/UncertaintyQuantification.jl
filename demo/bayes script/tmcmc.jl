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

function tresample(sset, hatw)
    snum = size(sset,1)
    rsset = zeros(snum, size(sset,2))
    hatwn = zeros(snum, 1)
    number = [1:snum...]

    for i in 1:snum
        resample = StatsBase.sample(number, ProbabilityWeights(vec(hatw)))
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

        sset, w = tresample(sset, hatw)
        hatw = w./sum(w)

        sig = calcsigma(sset, hatw, gamma)

        for i in 1:snum
            mvpropdist = MvNormal(zeros(size(sset, 2)), sig)
            propdistsample() = rand(mvpropdist)
            propdistpdf(x) = pdf(mvpropdist, x)
            likelihoodbeta(x) = likelihood(x)^(betaj1)
            sset[i,:] = MH(likelihoodbeta, prior, propdistsample, propdistpdf, sset[i,:], 2,1)
        end
        betaj = betaj1
    end

    return sset
end

function samplemean(
    dsample::Array
    )
    sumx = sum(dsample, dims=1)
    mcit = size(dsample,1)
    return sumx./mcit
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
sig = [1 0.8; 0.8 10]
mvnormal = MvNormal(mu, sig)
mvnormal2 = MvNormal([3, -5], sig)
likenormal(x) = pdf(mvnormal, [x[1], x[2]])+pdf(mvnormal2, [x[1], x[2]])

mpriorsample(n) = rand(Uniform(-40,40),(n, 2))

#tmsample = tmcmc(likenormal, prior, mpriorsample, 10000, gamma)

#display(histogram(msample[:, 1]))
#display(scatter(tmsample[:,1], tmsample[:, 2]))



function exlikelihood(x)
    eigvaln = [1.51 0.33; 4.01 0.30; 3.16 0.27; 3.21 0.18; 2.19 0.33; 1.71 0.23;
 2.73 0.21; 5.51 0.20; 1.95 0.11; 4.48 0.20; 1.43 0.16; 2.91 0.26; 3.81 0.23;
 3.58 0.25; 2.62 0.25]

    exsig = [1 0.5]
    rv = 0

    la1 = (x[1] + 2*x[2] + sqrt(x[1].^2 + 4x[2].^2))/2
    la2 = (x[1] + 2*x[2] - sqrt(x[1].^2 + 4x[2].^2))/2

    model = [la1 la2]


    for i in 1:2
        for j in 1:15
            rv += ((eigvaln[j,i] - model[i])/exsig[i])^2
        end
    end

    rv = exp(-1/2*rv)
    return rv
 end

 
 tuningsig = [0.04 0; 0 0.04]
 expropdist = MvNormal(mu, tuningsig)
 expropdistsample() = rand(expropdist)
 expropdistpdf(x) = pdf(expropdist, x)
 expriorsample(n) = rand(Uniform(0.001,4),(n, 2))

 function uniprior(x) 
    if x[1] > 0.001 && x[1] < 4 && x[2] > 0.001 && x[2] < 4
        rv = 1
    else 
        rv = 0
    end
    return  rv
end

 exsample = tmcmc(exlikelihood, uniprior, expriorsample, 5000, 0.1)
 
 display(scatter(exsample[:,1], exsample[:, 2]))
 #display(histogram(exsample[:, 2]))
