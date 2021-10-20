using LinearAlgebra, DataFrames, FiniteDifferences, Dierckx, Reexport,
Distributions, Random, Plots, StatsBase

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

#rsample = MH(expfunc, prior, propdistpdf, propdistsample, [3], mcit, 500)

#display(histogram(rsample[:, 1]))

#eqt = calcequaltails(rsample, 0.8)

mu = [12, 7]
sig = [4 0; 0 4]
mvnormal = MvNormal(mu, sig)
likenormal(x) = pdf(mvnormal, [x[1], x[2]])

propsig = [0.5 0; 0 0.5]
mpropdist = MvNormal([0;0], propsig)
mpropdistsample() = rand(mpropdist)
mpropdistpdf(x) = pdf(mpropdist, x)

msample = MH(likenormal, prior, mpropdistpdf, mpropdistsample, [-10; -7], 15000, 0)

display(plot(1:15000, msample[:,1]))
#display(histogram(msample[:, 1]))
#display(scatter(msample[:,1], msample[:, 2]))
#eqt = calcequaltails(msample, 0.8)



# Eigenvalue problem


 function exlikelihood(x)
    eigvaln = [1.51 0.33; 4.01 0.30; 3.16 0.27; 3.21 0.18; 2.19 0.33; 1.71 0.23;
 2.73 0.21; 5.51 0.20; 1.95 0.11; 4.48 0.20; 1.43 0.16; 2.91 0.26; 3.81 0.23;
 3.58 0.25; 2.62 0.25]

    exsig = [1 0.5]
    rv = 0

    la1 = (x[1] + 2*x[2] + sqrt(x[1].^2 + 4*x[2].^2))/2
    la2 = (x[1] + 2*x[2] - sqrt(x[1].^2 + 4*x[2].^2))/2

    model = [la1 la2]

    for i in 1:2
        for j in 1:15
            rv += ((eigvaln[j,i] - model[i])/exsig[i])^2
        end
    end

    return exp(-rv/2)
 end

 
 tuningsig = [0.04 0; 0 0.04]
 expropdist = MvNormal([0; 0], tuningsig)
 expropdistsample() = rand(expropdist)
 expropdistpdf(x) = pdf(expropdist, x)
 function uniprior(x) 
    if x[1] > 0.001 && x[1] < 4 && x[2] > 0.001 && x[2] < 4
        rv = 1
    else 
        rv = 0
    end
    return  rv
end

 #exsample = MH(exlikelihood, uniprior, expropdistpdf, expropdistsample, [2.84; 2.33], 5000, 0)
 
 #display(scatter(exsample[:,1], exsample[:, 2]))

 println(exlikelihood([0.15 1.5]))
 println(exlikelihood([3 2.8]))
 println(exlikelihood([0.5 1.5]))
 println(exlikelihood([1.3 1.0]))

 # change order sample pdf