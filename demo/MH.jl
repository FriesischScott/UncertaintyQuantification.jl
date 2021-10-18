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

    sample = zeros(mcit,size(startval,1))
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

function calcequaltails(
    dsample,
    quantile::Float64,
    )
    snum = size(dsample, 1)
    dsample = sort(dsample, dims = 1)

    low = floor(Int32,snum*(1-quantile)/2)
    high = floor(Int32,snum - snum*(1-quantile)/2)

    rv = zeros(2, size(dsample,2))

    for i in 1:size(dsample,2)
        rv[:,i] = [dsample[low,i]; dsample[high,i]]
    end
    return rv
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
propdist = Normal(0, 1)
propdistsample() = [rand(propdist)]
propdistpdf(x) = pdf(propdist, x[1])



mcit = 50500

#rsample = MH(expfunc, prior, propdistsample, propdistpdf, [3], mcit, 500)

#display(histogram(rsample[:, 1]))

#eqt = calcequaltails(rsample, 0.8)

mu = [0, 0]
sig = [400 0; 0 400]
mvnormal = MvNormal(mu, sig)
likenormal(x) = pdf(mvnormal, [x[1], x[2]])

propsig = [400 0; 0 400]
mpropdist = MvNormal(mu, propsig)
mpropdistsample() = rand(mpropdist)
mpropdistpdf(x) = pdf(mpropdist, x)

msample = MH(likenormal, prior, mpropdistsample, mpropdistpdf, [2; 4], 50500, 500)

display(histogram(msample[:, 1]))
#display(scatter(msample[:,1], msample[:, 2]))
#eqt = calcequaltails(msample, 0.8)
