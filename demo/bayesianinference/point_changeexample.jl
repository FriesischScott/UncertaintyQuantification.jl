using UncertaintyQuantification, Random

Random.seed!(2547);
N = 100
a = 2
b = 1

changepoint = floor(Int, N*rand(Uniform(0,1)))

lambda = Gamma(a, 1/b)

lambdas = zeros(N)

for i in 1:N
    if i<changepoint
        lambdas[i] = rand(lambda) *changepoint
    else
        lambdas[i] = rand(lambda) *(N-changepoint)
    end
end

sampledata = [rand(Poisson(i)) for i in lambdas]
Random.seed!()

x = sampledata


lambda1(in) = rand(Gamma(a + sum(x[1:changepoint-1]), 1/(changepoint + b)))
lambda2(in) = rand(Gamma(a + sum(x[changepoint:N]), 1/(N - changepoint + b)))

function multinomialN(in)
    propn = zeros(N)
    for i in 1:N
    propn[i] = sum(x[1:i-1])*log(in[1]) - i*in[1] + sum(x[i:N])*log(in[2]) - (N-changepoint)*in[2]
    end
    propn = exp.(propn .- maximum(propn))
    multiN = rand(Multinomial(1, propn./sum(propn)))
    rv = 0
    for i in 1:N
        if multiN[i]==1 rv = i end 
    end
    return rv
end

condist = [lambda1, lambda2, multinomialN]

rsample = gibbssample(condist, [3, 6, 50], 10000, 200)