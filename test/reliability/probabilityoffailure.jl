@testset "Probability of Failure" begin

    @testset "Monte Carlo" begin

        Random.seed!(8128)

        x = RandomVariable(Uniform(0.0, 1.0), :x)
        y = RandomVariable(Uniform(0.0, 1.0), :y)

        pf, _ = probability_of_failure(
            df -> 1 .- sqrt.(df.x.^2 + df.y.^2),
            [x, y],
            MonteCarlo(1000000),
        )

        Random.seed!()

        @test round(4.0 * (1 - pf), digits=3) == 3.141
    end

    @testset "Line sampling" begin

        Random.seed!(8128)

        x = RandomVariable(Uniform(0.0, 1.0), :x)
        y = RandomVariable(Uniform(0.0, 1.0), :y)
        r = Parameter(1, :r)

        d = Model(df -> sqrt.(df.x.^2 + df.y.^2), :d)

        g = df -> df.r .- df.d # performance function

        pf, _ = probability_of_failure(
            [d], g, [x, y, r], LineSampling(100),
        )

        Random.seed!()

        @test round(4.0 * (1 - pf), digits=3) == 3.142

        @test_logs (:warn, "All samples for line 1 are outside the failure domain") probability_of_failure(
            [d], g, [x, y, r], LineSampling(1, [0, 0.1]))
        @test_logs (:warn, "All samples for line 1 are inside the failure domain") probability_of_failure(
            [d], g, [x, y, r], LineSampling(1, [10, 20]))
    end

    @testset "Subset Simulation" begin
        # Kontantin Zuev - Subset Simulation Method for Rare Event Estimation: An Introduction
        # Example 6.1
        # Target pf of 1e-10

        Random.seed!(8128)

        pf_analytical = 1e-10

        x1 = RandomVariable(Normal(), :x1)
        x2 = RandomVariable(Normal(), :x2)

        y = Parameter(sqrt(2) * quantile(Normal(), 1 - pf_analytical), :y)

        g = Model(df -> df.x1 + df.x2, :g)

        F = df -> df.y .- df.g

        subset = SubSetSimulation(10^3, 0.1, 10, Normal())

        pf, _ = probability_of_failure(
            g,
            F,
            [x1, x2, y],
            subset
        )

        Random.seed!()

        @test isapprox(pf, pf_analytical, atol=10^-10)
    end
end