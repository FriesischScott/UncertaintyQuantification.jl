using UncertaintyQuantification

function likelihood(x)
    eigvaln = [1.51 0.33; 4.01 0.30; 3.16 0.27; 3.21 0.18; 2.19 0.33; 1.71 0.23;
 2.73 0.21; 5.51 0.20; 1.95 0.11; 4.48 0.20; 1.43 0.16; 2.91 0.26; 3.81 0.23;
 3.58 0.25; 2.62 0.25]

    sig = [1 0.5]

    la1 = (x[1] + 2*x[2] + sqrt(x[1].^2 + 4x[2].^2))/2
    la2 = (x[1] + 2*x[2] - sqrt(x[1].^2 + 4x[2].^2))/2

    model = [la1 la2]
    rv = 0

    for i in 1:2
        for j in 1:15
            rv += ((eigvaln[j,i] - model[i])/sig[i])^2
        end
    end

    rv = exp(-1/2*rv)
    return rv
 end

 
 tuningsig = [0.04 0; 0 0.04]
 propdist = MvNormal([0;0], tuningsig)
 propdistsample() = rand(expropdist)
 propdistpdf(x) = pdf(expropdist, x)

 priorsample(n) = rand(Uniform(0.001,4),(n, 2))

 function uniformprior(x) 
    if x[1] > 0.001 && x[1] < 4 && x[2] > 0.001 && x[2] < 4
        rv = 1
    else 
        rv = 0
    end
    return  rv
end



tmcmcsample = tmcmc(likelihood, uniformprior, priorsample, 10000, 0.1)