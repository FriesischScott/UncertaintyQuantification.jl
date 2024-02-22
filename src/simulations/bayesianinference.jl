function mh(
    likelihood::Function,
    prior::Function,
    propdistpdf::Function,
    propdistsample::Function,
    startval::Vector,
    mcit::Int64,
    burn::Int64
    )

    sample = zeros(mcit,length(startval))
    sample[1,:] = startval

    for i in 1:mcit-1
        curr = sample[i,:]

        #create new sample
        proposalval = curr + propdistsample()

        #calculate acceptance probability
        A = (likelihood(proposalval)/likelihood(curr)*prior(proposalval)/
            prior(curr)*propdistpdf(proposalval-curr)/
            propdistpdf(curr-proposalval))[1]

        #accept or reject sample
        if A>=rand(Uniform(0,1))
            sample[i+1,:] = proposalval
        else
            sample[i+1,:] = curr
        end
    end
    return sample[burn+1:mcit,:]
end

function gibbssample(
    condist::Vector,
    startval::Vector,
    mcit::Int64,
    burn::Int64
    )

    sample = zeros(mcit, length(startval))
    sample[1,:] = startval

    for i in 1:mcit-1
        sample[i+1,:] = sample[i,:]
        for j in 1:size(sample, 2)
            #update each variable by itself
            sample[i+1,j] = condist[j](sample[i+1,:])
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
        #resample according to weight
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

    betaj1 = betaj + (1.1 - betaj)*0.5
    likelihoodcov = zeros(snum, 1)
    step = (1.1 - betaj)*0.25
    cov = 0

    for i in 1:snum
        likelihoodcov[i] = likelihood(sset[i,:])
    end
    #beta approximates so that coefficient of variance --> 1
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
    gamma::Real
    )
    sdim = size(sset, 2)
    snum = size(sset, 1)

    thet = zeros(sdim, 1)
    sig = zeros(sdim, sdim)

    #formula for calculating Sigma in Ching and Chen 2007
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
    gamma::Real
    )

    j = 0
    betaj = 0.0
    sset = priorsample(snum)


    while betaj < 1
        betaj1, w = calcbeta(likelihood, betaj, sset, snum)
        hatw = w./sum(w)#normalize

        sset, w = tresample(sset, hatw)
        hatw = w./sum(w)

        sig = calcsigma(sset, hatw, gamma)

        for i in 1:snum
            #MH step for each sample
            mvpropdist = MvNormal(zeros(size(sset, 2)), sig)
            propdistsample() = rand(mvpropdist)
            propdistpdf(x) = pdf(mvpropdist, x)
            likelihoodbeta(x) = likelihood(x)^(betaj1)
            sset[i,:] = mh(likelihoodbeta, prior, propdistpdf, propdistsample, sset[i,:], 2,1)
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
        #weighted resampling
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
    term::Function
    )

    j = 0
    sset = priorsample(snum)
    w = zeros(snum,1)

    #calculate weights
    for i in 1:snum
        w[i] = likelihood(sset[i,:])*prior(sset[i,:])
    end

    hatw = w./sum(w)#normalize

    while term(sset, j)
        j += 1
        Neff = 1/sum(hatw.^2)#efficient sample size
        if Neff < snum/2 #test for degeneracy
            sset = resample(sset, hatw)
            hatw = ones(snum,1)./snum
        end

        rsset = zeros(snum, size(sset,2))
        w = zeros(snum, 1)

        for i in 1:snum
            #one step in a random direction
            rsset[i,:] = sset[i,:] + propdistsample()
            w[i] = hatw[i].*likelihood(rsset[i,:])/likelihood(sset[i,:])*
                prior(rsset[i,:])/prior(sset[i,:])*
                propdistpdf(rsset[i,:]-sset[i,:])/propdistpdf(sset[i,:]-rsset[i,:])
        end

        sset = rsset
        hatw = w./sum(w)
    end
    return sset
end

function calcequaltails(
    dsample::Array,
    quantile::Float64,
    )
    snum = size(dsample, 1)
    dsample = sort(dsample, dims = 1)

    #calc cutoff points
    low = floor(Int32,snum*(1-quantile)/2)
    high = floor(Int32,snum - snum*(1-quantile)/2)

    rv = zeros(2, size(dsample,2))

    for i in 1:size(dsample,2)
        rv[:,i] = [dsample[low,i]; dsample[high,i]]
    end
    return rv
end

function grconvergence(
    likelihood::Function,
    prior::Function,
    priorsample::Function,
    propdistpdf::Function,
    propdistsample::Function,
    mcit::Int64,
    nchains::Int64
    )

    startval = priorsample(nchains)
    dims = size(startval, 2)
    rhat = 0
    mcitb = floor(Int, mcit/2)
    samples = zeros(floor(Int, mcit/2), size(startval,2), nchains)

    #create n Markov chains
    for i in 1:nchains
        samples[:,:,i] = mh(likelihood, prior, propdistpdf, propdistsample,startval[i,:], mcit, floor(Int, mcit/2))
    end

    if dims == 1
        chainvariance = zeros(nchains)
        chainmean = zeros(nchains)
        overallmean = 0
        for i in 1:nchains
            chainmean[i] = mean(samples[:,1,i], dims=1)[1]
            overallmean += chainmean[i]/nchains
            for j in 1:floor(Int,mcit/2) chainvariance[i] += (samples[j,1,i] - chainmean[i])^2 end
        end
        #vector of variances for each chain
        chainvariance = chainvariance./(mcitb-1)

        #variance of all chains
        mvariance = 0
        for i in 1:nchains mvariance += chainvariance[i]/nchains end

        #variance bias estimator
        Bn = 0
        for i in 1:nchains Bn += (chainmean[i]-overallmean)^2 end
        Bn = Bn/(nchains-1)

        #over corrected variance
        sig2 = ((mcitb-1)/mcitb)*mvariance + Bn

        rhat = sqrt(sig2/mvariance)
    else
        chainvariance = zeros(dims, dims, nchains)
        chainmean = zeros(dims, nchains)
        overallmean = zeros(dims)

        for i in 1:nchains
            chainmean[:, i] = mean(samples[:,:,i], dims=1)
            overallmean += chainmean[:, i]./nchains
            for j in 1:floor(Int, mcit/2)
                v = samples[j,:,i]-chainmean[:,i]
                chainvariance[:,:,i] += v*v'
            end
        end
        #array of covariance matrices for each chain
        chainvariance = chainvariance./(mcitb-1)

        #covariance matrix of all chains
        mvariance = zeros(dims, dims)
        for i in 1:nchains mvariance += chainvariance[:,:,i]./nchains end

        #bias estimator
        Bn = zeros(dims, dims)
        for i in 1:nchains
            v = chainmean[:,i] - overallmean
            Bn += v*v'
        end

        Bn = Bn*mcitb/(nchains-1)
        mvariance = inv(mvariance)
        maxeig = maximum(eigvals(mvariance*Bn))

        rhat = sqrt((mcitb-1)/mcitb + maxeig/mcitb)
    end

    return rhat
end
