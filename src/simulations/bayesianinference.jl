function MH(
    likelihood::Function,
    prior::Function,
    propdistpdf::Function,
    propdistsample::Function,
    startval::Array,
    mcit::Int64,
    burn::Int64
    )

    sample = zeros(mcit,size(startval,1))
    sample[1,:] = startval

    for i in 1:mcit-1
        curr = sample[i,:]

        proposalval = curr + propdistsample()

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

function gibbssample(
    condist::Array,
    startval::Array,
    mcit::Int64,
    burn::Int64
    )

    sample = zeros(mcit, size(startval, 2))
    sample[1,:] = startval

    for i in 1:mcit-1
        sample[i+1,:] = sample[i,:]
        for j in 1:size(sample, 2)
            sample[i+1,j] = rand(condist[j](sample[i+1,:]))
        end
    end
    return sample[burn+1:mcit,:]
end

function tresample(
    sset::Array,
    hatw::Array
    )
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
    sset::Array,
    hatw::Array,
    gamma::Int64
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
            sset[i,:] = MH(likelihoodbeta, prior, propdistpdf, propdistsample, sset[i,:], 2,1)
        end
        betaj = betaj1
    end

    return sset
end

function resample(
    sset::Array,
    hatw::Array
    )
    snum = size(sset,1)
    rsset = zeros(snum, size(sset,2))
    number = [1:snum...]

    for i in 1:snum
        resample = StatsBase.sample(number, ProbabilityWeights(vec(hatw)))
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
