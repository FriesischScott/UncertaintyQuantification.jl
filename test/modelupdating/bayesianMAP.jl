
@testset "ML and MAP estimates" begin
    N_binom = 15
    data_binom = [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1] # p = 0.8

    N_normal = 25
    data_normal_var = [
        3.1170072530410713,
        5.793425513286137,
        3.019644523423345,
        4.198492885643784,
        4.3404411932408005,
        5.338786752835615,
        1.6563692415479587,
        3.496696575234516,
        5.093428444721792,
        5.091321607109399,
        4.994437068835541,
        4.285932709130585,
        1.6751627811539134,
        5.51537861306012,
        2.4998958266295026,
        7.065078073070958,
        5.005937831600564,
        3.9349403146434914,
        3.9340063190556913,
        4.400921498325048,
        5.276032641042012,
        4.034478885301326,
        7.5811254216258215,
        6.468279762464751,
        4.416669261547403,
    ] # μ = 4, σ² = 5

    data_normal_mean = [
        3.2447586269995603,
        7.31570639905421,
        2.9583720841640524,
        2.7006645805494793,
        -5.390670487961957,
        13.222284857030616,
        4.274588755848634,
        6.92516681771173,
        -14.259515860959095,
        2.8938242850671676,
        8.702867817886007,
        17.60466794732554,
        -3.791069613048381,
        3.8733459676840587,
        5.010069312798485,
        -1.8534537575475099,
        10.843101672796852,
        -1.7035437372547804,
        -1.8363450654785087,
        -2.929811781116549,
        -1.6820361524479646,
        -1.5074002658326515,
        1.2143899492909962,
        2.9977426550531416,
        3.9469743716484778,
    ] # μ = 4, σ = 5

    function binomialinferencebenchmark(
        estimater::AbstractBayesianPointEstimate, prior::Beta{Float64}
    )
        alpha0 = prior.α
        beta0 = prior.β

        alpha_posterior = alpha0 + sum(data_binom)
        beta_posterior = beta0 + N_binom - sum(data_binom)

        analytic_mode = (alpha_posterior-1) / (alpha_posterior + beta_posterior - 2)
        analytic_mean = alpha_posterior / (alpha_posterior + beta_posterior)
        analytic_var =
            alpha_posterior * beta_posterior /
            ((alpha_posterior + beta_posterior)^2 * (alpha_posterior + beta_posterior + 1))

        function loglikelihood(df)
            return [
                sum(logpdf.(Binomial.(N_binom, df_i.x), sum(data_binom))) for
                df_i in eachrow(df)
            ]
        end

        logprior(df) = logpdf.(prior, df.x)

        estimate = bayesianupdating(loglikelihood, UQModel[], estimater, prior = logprior)

        return estimate, analytic_mode, analytic_mean
    end

    function normalmeanbenchmark(estimater::AbstractBayesianPointEstimate, prior::Normal{Float64})
        std_fixed = 5     # Fixed

        prior_mean = prior.μ
        prior_std = prior.σ

        analytic_std = 1 / ((1 / prior_std^2) + (N_normal / std_fixed^2))
        analytic_mean =
            analytic_std * (prior_mean / prior_std^2 + sum(data_normal_mean) / std_fixed^2)

        function loglikelihood(df)
            return [
                sum(logpdf.(Normal.(df_i.x, std_fixed), data_normal_mean)) for
                df_i in eachrow(df)
            ]
        end

        logprior(df) = logpdf.(prior, df.x)

        estimation = bayesianupdating(loglikelihood, UQModel[], estimater, prior = logprior)

        return estimation, analytic_mean, analytic_std
    end

    @testset "MAP normal mean analytical" begin

        optMethod = "LBFGS"
        x0 = [0.]

        prior_mean = 2
        prior_std = 10

        prior = Normal(prior_mean, prior_std)

        prior_sample_ = RandomVariable(prior, :x)

        estimater = MaximumAPosterioriBayesian([prior_sample_], optMethod, x0)

        MAPEst, analytic_mean, analytic_std = normalmeanbenchmark(estimater, prior)

        @test mean(MAPEst.x) ≈ analytic_mean rtol = 0.05

    end

    @testset "MLE normal mean analytical" begin

        optMethod = "LBFGS"
        x0 = [0.]

        prior_mean = 2
        prior_std = 10

        prior = Normal(prior_mean, prior_std)

        prior_sample_ = RandomVariable(prior, :x)

        estimater = MaximumLikelihoodBayesian([prior_sample_], optMethod, x0)

        MLEst, analytic_mean, analytic_std = normalmeanbenchmark(estimater, prior)

        @test mean(MLEst.x) ≈ analytic_mean rtol = 0.05

    end

    @testset "MAP binomial test" begin

        optMethod = "LBFGS"
        x0 = [0.5]
        
        prior_Function = Beta(1,1)
        prior = RandomVariable(prior_Function, :x)

        estimater = MaximumAPosterioriBayesian([prior], optMethod, x0; lowerbounds = [0.], upperbounds = [1.])

        estimate, analytic_mode, analytic_mean = binomialinferencebenchmark(estimater, prior_Function)

        @test mean(estimate.x) ≈ analytic_mode rtol = 0.1
    end

    @testset "MLE binomial test" begin

        optMethod = "LBFGS"
        x0 = [0.5]
        
        prior_Function = Beta(1,1)
        prior = RandomVariable(prior_Function, :x)

        estimater = MaximumLikelihoodBayesian([prior], optMethod, x0; lowerbounds = [0.], upperbounds = [1.])

        estimate, analytic_mode, analytic_mean = binomialinferencebenchmark(estimater, prior_Function)

        @test mean(estimate.x) ≈ analytic_mode rtol = 0.1
    end

end
