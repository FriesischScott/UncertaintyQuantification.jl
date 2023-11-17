@testset "Probability of Failure" begin
    @testset "Monte Carlo" begin
        x = RandomVariable(Uniform(0.0, 1.0), :x)
        y = RandomVariable(Uniform(0.0, 1.0), :y)

        pf, _ = probability_of_failure(
            df -> 1 .- sqrt.(df.x .^ 2 + df.y .^ 2), [x, y], MonteCarlo(1000000)
        )

        pi = 4.0 * (1 - pf)

        @test pi ≈ 3.141 atol = 0.01
    end

    @testset "Importance Sampling" begin
        # Englund and Rackwitz - A benchmark study on importance sampling techniques
        # in structural reliability (1993)
        # Example 1
        u = RandomVariable.(Normal(), [:u1, :u2])

        g = df -> 2^(1 / 2) .- (df.u1 .+ df.u2)

        _, β, dp, α = probability_of_failure(g, u, FORM())

        pf, _, _ = probability_of_failure(g, u, ImportanceSampling(10000, β, dp, α))

        @test pf ≈ 0.159 rtol = 0.05
    end

    @testset "Line sampling" begin
        # Englund and Rackwitz - A benchmark study on importance sampling techniques
        # in structural reliability (1993)
        # Example 1
        u = RandomVariable.(Normal(), [:u1, :u2])

        g = df -> 2^(1 / 2) .- (df.u1 .+ df.u2)

        pf, _ = probability_of_failure(g, u, LineSampling(100))

        @test pf ≈ 0.159 rtol = 0.05

        @test_logs (:warn, "All samples for line 1 are outside the failure domain") probability_of_failure(
            g, u, LineSampling(1, [0, 0.1])
        )
        @test_logs (:warn, "All samples for line 1 are inside the failure domain") probability_of_failure(
            g, u, LineSampling(1, [10, 20])
        )
    end

    # Kontantin Zuev - Subset Simulation Method for Rare Event Estimation: An Introduction
    # Example 6.1
    # Target pf of 1e-10
    pf_analytical = 1e-10

    x1 = RandomVariable(Normal(), :x1)
    x2 = RandomVariable(Normal(), :x2)

    y = Parameter(sqrt(2) * quantile(Normal(), 1 - pf_analytical), :y)

    g = Model(df -> df.x1 + df.x2, :g)

    F = df -> df.y .- df.g

    @testset "Subset Simulation" begin
        subset = SubSetSimulation(10^4, 0.1, 20, Normal())

        pf, _, _ = probability_of_failure(g, F, [x1, x2, y], subset)

        # 95% conf intervals estimated from 1000 runs
        @test 1.4e-11 < pf < 9.1e-10
    end

    @testset "Subset Infinity" begin
        subset = SubSetInfinity(10^4, 0.1, 20, 0.5)

        pf, _, _ = probability_of_failure(g, F, [x1, x2, y], subset)

        # 95% conf intervals estimated from 1000 runs
        @test 1.8e-11 < pf < 1.6e-9
    end

    @testset "Subset Infinity Adaptive" begin
        subset = SubSetInfinityAdaptive(10^4, 0.1, 20, 100, 0.6, 1.0)

        pf, _, _ = probability_of_failure(g, F, [x1, x2, y], subset)

        # 95% conf intervals estimated from 1000 runs
        @test 3.14e-11 < pf < 3.4e-10
    end
end
