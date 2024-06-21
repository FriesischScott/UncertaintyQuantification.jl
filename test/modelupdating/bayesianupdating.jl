@testset "Bayesian Updating" begin
    @testset "Single Component Metropolis Hastings" begin
        prior(df) = pdf.(Normal(), df.x)
        likelihood(df) = pdf.(Normal(2, 0.5), df.x)

        proposal = Normal()
        x0 = (x=0.0,)
        n = 1000
        burnin = 100

        mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)
        mh_samples, _ = bayesianupdating(prior, likelihood, mh)

        @test mean(mh_samples.x) ≈ 1.6 rtol = 0.05
        @test std(mh_samples.x) ≈ 0.45 rtol = 0.05
    end

    @testset "Transitional Markov Chain Monte Carlo" begin
        prior_sample = [RandomVariable(Normal(), :x)]
        prior(df) = logpdf.(Normal(), df.x)

        likelihood(df) = logpdf.(Normal(2, 0.5), df.x)

        n = 1000
        burnin = 100

        tmcmc = TransitionalMarkovChainMonteCarlo(prior_sample, n, burnin, 0.2)

        tmcmc_samples, _ = bayesianupdating(prior, likelihood, UQModel[], tmcmc)

        @test mean(tmcmc_samples.x) ≈ 1.6 rtol = 0.05
        @test std(tmcmc_samples.x) ≈ 0.45 rtol = 0.05
    end

    @testset "TMCMC binomal inference analytical" begin
        p_true = 0.8
        N_data = 15

        Data = rand(N_data) .<= p_true

        alpha0 = 1
        beta0 = 1

        # Posterior beta parameters
        alpha_posterior = alpha0 + sum(Data)
        beta_posterior = beta0 + N_data - sum(Data)
        # posterior_exact(x) = pdf(Beta(alpha_posterior, beta_posterior), x)

        analytic_mean = alpha_posterior/(alpha_posterior + beta_posterior)
        analytic_var = alpha_posterior*beta_posterior/((alpha_posterior+beta_posterior)^2 * (alpha_posterior+beta_posterior+1))

        logprior(df) = logpdf(Beta(alpha0, beta0), df.x)
        loglikelihood(df) = [sum(logpdf.(Binomial.(N_data, df_i.x), sum(Data))) for df_i in eachrow(df)]
        prior_sample_ = RandomVariable(Beta(alpha0, beta0), :x)

        n = 2000
        burnin = 200
        thin = 5

        tmcmc = TransitionalMarkovChainMonteCarlo([prior_sample_], n, burnin, thin)
        tmcmc_samples, _ = bayesianupdating(logprior, loglikelihood, UQModel[], tmcmc)

        @test mean(tmcmc_samples.x) ≈ analytic_mean rtol = 0.05
        @test var(tmcmc_samples.x) ≈ analytic_var rtol = 0.05
    end

    @testset "TMCMC normal mean inference analytical" begin
        mean_true = 4
        std_fixed = 5     # Fixed

        N_data = 25

        Data = rand(Normal(mean_true, std_fixed), N_data)

        prior_mean = 2
        prior_std = 10

        # Exact conjugate posterior moments
        analytic_std = 1/((1/prior_std^2) + (N_data/std_fixed^2))
        analytic_mean = analytic_std * (prior_mean / prior_std^2 + sum(Data) /std_fixed^2)
        # posterior_exact(x) = pdf(Normal(analytic_mean, analytic_std), x)

        prior_sample = [RandomVariable(Normal(), :x)]
        prior(df) = logpdf.(Normal(), df.x)

        likelihood(df) = logpdf.(Normal(2, 0.5), df.x)

        logprior(df) = logpdf(Normal(prior_mean, prior_std), df.x)
        prior_sample_ = RandomVariable(Normal(prior_mean, prior_std), :x)
        loglikelihood(df) = [sum(logpdf.(Normal.(df_i.x, std_fixed), Data)) for df_i in eachrow(df)]
        
        
        n = 2000
        burnin = 200
        thin = 5
        
        tmcmc = TransitionalMarkovChainMonteCarlo([prior_sample_], n, burnin, thin)
        tmcmc_samples, _ = bayesianupdating(logprior, loglikelihood, UQModel[], tmcmc)
        
        @test mean(tmcmc_samples.x) ≈ analytic_mean rtol = 0.05
        @test std(tmcmc_samples.x) ≈ analytic_std rtol = 0.05
    end
end
