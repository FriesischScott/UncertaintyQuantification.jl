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

    @testset "Line sampling" begin
        x = RandomVariable(Uniform(0.0, 1.0), :x)
        y = RandomVariable(Uniform(0.0, 1.0), :y)
        r = Parameter(1, :r)

        d = Model(df -> sqrt.(df.x .^ 2 + df.y .^ 2), :d)

        g = df -> df.r .- df.d # performance function

        pf, _ = probability_of_failure([d], g, [x, y, r], LineSampling(100))

        @test 4.0 * (1 - pf) ≈ 3.141 atol = 0.01

        @test_logs (:warn, "All samples for line 1 are outside the failure domain") probability_of_failure(
            [d], g, [x, y, r], LineSampling(1, [0, 0.1])
        )
        @test_logs (:warn, "All samples for line 1 are inside the failure domain") probability_of_failure(
            [d], g, [x, y, r], LineSampling(1, [10, 20])
        )
    end

    @testset "Subset Simulation" begin
        # Kontantin Zuev - Subset Simulation Method for Rare Event Estimation: An Introduction
        # Example 6.1
        # Target pf of 1e-10
        pf_analytical = 1e-10

        x1 = RandomVariable(Normal(), :x1)
        x2 = RandomVariable(Normal(), :x2)

        y = Parameter(sqrt(2) * quantile(Normal(), 1 - pf_analytical), :y)

        g = Model(df -> df.x1 + df.x2, :g)

        F = df -> df.y .- df.g

        subset = SubSetSimulation(10^4, 0.1, 20, Normal())

        pf, _, _ = probability_of_failure(g, F, [x1, x2, y], subset)

        # 95% conf intervals estimated from 1000 runs
        @test 1.4e-11 < pf < 9.1e-10

        subset = SubSetInfinity(10^4, 0.1, 20, 0.5)

        pf, _, _ = probability_of_failure(g, F, [x1, x2, y], subset)

        # 95% conf intervals estimated from 1000 runs
        @test 1.8e-11 < pf < 1.6e-9
    end
end
