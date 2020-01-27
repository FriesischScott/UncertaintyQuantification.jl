@testset "Probability of Failure" begin

    @testset "Monte Carlo" begin

        Random.seed!(8128)

        x = RandomVariable(Uniform(0.0, 1.0), :x)
        y = RandomVariable(Uniform(0.0, 1.0), :y)

        d = Model(df -> sqrt.(df.x .^ 2 + df.y .^ 2), :d)

        pf, _ = probability_of_failure(
            [d],
            df -> 1 .- df.d,
            [x, y],
            MonteCarlo(1000000),
        )

        Random.seed!()

        @test round(4.0 * (1 - pf), digits = 3) == 3.141
    end

    @testset "Line sampling" begin

        Random.seed!(8128)

        x = RandomVariable(Uniform(0.0, 1.0), :x)
        y = RandomVariable(Uniform(0.0, 1.0), :y)

        d = Model(df -> sqrt.(df.x .^ 2 + df.y .^ 2), :d)

        pf, _ = probability_of_failure(
            [d],
            df -> 1 .- df.d,
            [x, y],
            LineSampling(100),
        )

        Random.seed!()

        @test round(4.0 * (1 - pf), digits = 3) == 3.142
    end

end
