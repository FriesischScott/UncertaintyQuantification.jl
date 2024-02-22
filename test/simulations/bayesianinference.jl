sig = [0.5 0; 0 0.5]
mvnormal = MvNormal([3, 5], sig)
mvnormal2 = MvNormal([2, 3], sig)
likelihood(x) = pdf(mvnormal, [x[1], x[2]]) + pdf(mvnormal2, [x[1], x[2]])

prior(x) = 1
priorsample(n) = rand(Uniform(0.001,10),(n, 2))

propsig = [0.2 0; 0 0.2]
propdist = MvNormal([0,0], propsig)
propdistsample() = rand(propdist)
propdistpdf(x) = pdf(propdist, x)

function term(sset, j)
    rv = true
    if j >= 1 rv = false end
    return rv
 end

@testset "Metropolis-Hasting" begin
    Random.seed!(8128)
    mhsample = mh(likelihood, prior, propdistpdf, propdistsample, [3, 1], 15000, 500)
    Random.seed!()
    
    @test size(mhsample, 1) == 14500
    @test isapprox(mean(mhsample, dims=1), [2.5 4], atol=0.05)
end

@testset "TMCMC" begin
    Random.seed!(8128)
    tmcmcsample = tmcmc(likelihood, prior, priorsample, 15000, 0.2)
    Random.seed!()

    @test size(tmcmcsample, 1) == 15000
    @test isapprox(mean(tmcmcsample, dims=1), [2.5 4], atol=0.05)
end

@testset "SMC" begin
    Random.seed!(8128)
    sample = smc(likelihood, prior, priorsample, propdistpdf, propdistsample, 15000, term)
    Random.seed!()
    
    @test size(sample, 1) == 15000
    @test isapprox(mean(sample, dims=1), [2.5 4], atol=0.05)
end

@testset "Gelman-Rubin_convergence" begin
    Random.seed!(8128)
    rhat = grconvergence(likelihood, prior, priorsample, propdistpdf, propdistsample, 15000, 8)
    Random.seed!()
    
    @test isapprox(rhat, 1, atol=0.01)
end