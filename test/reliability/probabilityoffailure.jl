@testset "Probability of Failure" begin

    @testset "Monte Carlo" begin

        Random.seed!(8128)

        x = RandomVariable(Uniform(0.0, 1.0), :x)
        y = RandomVariable(Uniform(0.0, 1.0), :y)

        d = Model(df -> sqrt.(df.x.^2 + df.y.^2), :d)

        pf, _ = probability_of_failure(
            [d],
            df -> 1 .- df.d,
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

end
