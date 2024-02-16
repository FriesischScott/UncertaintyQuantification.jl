@testset "Bayesian Updating" begin
    @testset "Single Component Metropolis Hastings" begin
        prior(df) = pdf(Normal(), df.x)
        likelihood(df) = pdf(Normal(2, 0.5), df.x)

        proposal = Normal()
        x0 = (x=0.0,)
        n = 1000
        burnin = 100

        mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)
        mh_samples, _ = bayesianupdating(prior, likelihood, mh)

        @test mean(mh_samples.x) ≈ 1.6 rtol = 0.05
        @test std(mh_samples.x) ≈ 0.45 rtol = 0.05
    end
end
