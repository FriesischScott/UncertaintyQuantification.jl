using UncertaintyQuantification, Random
Random.seed!(2547);
N = 200
mu = [0, 0]
sig = [1 -0.55; -0.55 1]
mvnormal = MvNormal(mu, sig)


sam = zeros(200,2)
for i in 1:200 sam[i,:] = rand(mvnormal) end

Random.seed!()

function likelihood(x) 
    if x[1] >1 || x[1]<-1
        rv = 0
    else
        rv = 1
        for i in 1:N
            rv = rv/(2*pi*sqrt(1-x[1]^2))*exp(-1/(2*(1-x[1]^2))*(sam[i,1]^2-2*x[1]*
            sam[i,1]*sam[i,2]+sam[i,2]^2))
        end
    end
    return rv
end

function prior(x) 
    if x[1] >1 || x[1]<-1
        rv = 0
    else
        rv = 1/(1-x[1]^2)^(3/2)
    end
    return rv
end
priorsample(n) = rand(Uniform(-1,1),(n))


tmcmcsample = tmcmc(likelihood, prior, priorsample, 2000, 0.2)