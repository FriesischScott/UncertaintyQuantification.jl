@testset "Bayesian Updating" begin
    function binomialinferencebenchmark(
        sampler::AbstractBayesianMethod, prior::Beta{Float64}, p_true=0.8, N_data=15
    )
        alpha0 = prior.α
        beta0 = prior.β
        Data = rand(N_data) .<= p_true

        alpha_posterior = alpha0 + sum(Data)
        beta_posterior = beta0 + N_data - sum(Data)

        analytic_mean = alpha_posterior / (alpha_posterior + beta_posterior)
        analytic_var =
            alpha_posterior * beta_posterior /
            ((alpha_posterior + beta_posterior)^2 * (alpha_posterior + beta_posterior + 1))

        function loglikelihood(df)
            return [
                sum(logpdf.(Binomial.(N_data, df_i.x), sum(Data))) for df_i in eachrow(df)
            ]
        end

        logprior(df) = logpdf(prior, df.x)

        mcmc_samples, _ = bayesianupdating(logprior, loglikelihood, UQModel[], sampler)

        return mcmc_samples, analytic_mean, sqrt(analytic_var)
    end

    function normalmeanbenchmark(
        sampler::AbstractBayesianMethod, prior::Normal{Float64}, mean_true=4, N_data=25
    )
        std_fixed = 5     # Fixed

        Data = rand(Normal(mean_true, std_fixed), N_data)

        prior_mean = prior.μ
        prior_std = prior.σ

        analytic_std = 1 / ((1 / prior_std^2) + (N_data / std_fixed^2))
        analytic_mean = analytic_std * (prior_mean / prior_std^2 + sum(Data) / std_fixed^2)

        function loglikelihood(df)
            return [sum(logpdf.(Normal.(df_i.x, std_fixed), Data)) for df_i in eachrow(df)]
        end

        logprior(df) = logpdf(prior, df.x)

        mcmc_samples, _ = bayesianupdating(logprior, loglikelihood, UQModel[], sampler)

        return mcmc_samples, analytic_mean, analytic_std
    end

    function normalvarbenchmark(
        sampler::AbstractBayesianMethod, prior::InverseGamma, var_true=5, N_data=25
    )
        mean_fixed = 4     # Fixed

        Data = rand(Normal(mean_fixed, sqrt(var_true)), N_data)

        prior_shape = prior.invd.α
        prior_scale = prior.θ

        posterior_shape = prior_shape + N_data / 2
        posterior_scale = prior_scale + sum((Data .- mean_fixed) .^ 2) / 2

        posterior_exact = InverseGamma(posterior_shape, posterior_scale)

        function loglikelihood(df)
            return [
                sum(logpdf.(Normal.(mean_fixed, sqrt(df_i.x)), Data)) for
                df_i in eachrow(df)
            ]
        end

        function logprior1(df)  # Required because inverse gamma throws error for negative values
            if df.x < 0
                return -Inf
            end
            return logpdf(prior, df.x)
        end

        logprior(df) = [logprior1(df_i) for df_i in eachrow(df)]

        mcmc_samples, _ = bayesianupdating(logprior, loglikelihood, UQModel[], sampler)

        return mcmc_samples, mean(posterior_exact), std(posterior_exact)
    end

    @testset "Single Component MH binomal inference analytical" begin
        proposal = Normal()
        x0 = (x=0.5,)
        n = 2000
        burnin = 500

        prior = Beta(1, 1)

        mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

        mc_samples, analytic_mean, analytic_std = binomialinferencebenchmark(mh, prior)

        @test mean(mc_samples.x) ≈ analytic_mean rtol = 0.1
        @test std(mc_samples.x) ≈ analytic_std rtol = 0.1
    end

    @testset "Single Component MH normal mean analytical" begin
        proposal = Normal()
        x0 = (x=1.0,)
        n = 2000
        burnin = 500

        prior_mean = 2
        prior_std = 10

        prior = Normal(prior_mean, prior_std)

        mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

        mc_samples, analytic_mean, analytic_std = normalmeanbenchmark(mh, prior)

        @test mean(mc_samples.x) ≈ analytic_mean rtol = 0.1
        @test std(mc_samples.x) ≈ analytic_std rtol = 0.1
    end

    @testset "Single Component MH normal var analytical" begin
        proposal = Normal()
        x0 = (x=4.0,)
        n = 2000
        burnin = 500

        prior_shape = 30
        prior_scale = 100

        prior = InverseGamma(prior_shape, prior_scale)

        mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

        mc_samples, analytic_mean, analytic_std = normalvarbenchmark(mh, prior)

        @test mean(mc_samples.x) ≈ analytic_mean rtol = 0.1
        @test std(mc_samples.x) ≈ analytic_std rtol = 0.1
    end

    @testset "TMCMC binomal inference analytical" begin
        n = 2000
        burnin = 200

        prior = Beta(1, 1)

        prior_sample_ = RandomVariable(prior, :x)

        tmcmc = TransitionalMarkovChainMonteCarlo([prior_sample_], n, burnin)

        mc_samples, analytic_mean, analytic_std = binomialinferencebenchmark(tmcmc, prior)

        @test mean(mc_samples.x) ≈ analytic_mean rtol = 0.05
        @test std(mc_samples.x) ≈ analytic_std rtol = 0.05
    end

    @testset "TMCMC normal mean analytical" begin
        n = 2000
        burnin = 200

        prior_mean = 2
        prior_std = 10

        prior = Normal(prior_mean, prior_std)

        prior_sample_ = RandomVariable(prior, :x)

        tmcmc = TransitionalMarkovChainMonteCarlo([prior_sample_], n, burnin)

        mc_samples, analytic_mean, analytic_std = normalmeanbenchmark(tmcmc, prior)

        @test mean(mc_samples.x) ≈ analytic_mean rtol = 0.05
        @test std(mc_samples.x) ≈ analytic_std rtol = 0.05
    end

    @testset "TMCMC normal var analytical" begin
        n = 2000
        burnin = 200

        prior_shape = 30
        prior_scale = 100

        prior = InverseGamma(prior_shape, prior_scale)

        prior_sample_ = RandomVariable(prior, :x)

        tmcmc = TransitionalMarkovChainMonteCarlo([prior_sample_], n, burnin)

        mc_samples, analytic_mean, analytic_std = normalvarbenchmark(tmcmc, prior)

        @test mean(mc_samples.x) ≈ analytic_mean rtol = 0.05
        @test std(mc_samples.x) ≈ analytic_std rtol = 0.05
    end
end
