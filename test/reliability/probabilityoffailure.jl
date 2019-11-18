@testset "Probability of Failure" begin

    @testset "Monte Carlo" begin

        Random.seed!(8128)

        x = RandomVariable(Uniform(0, 1), :x)
        y = RandomVariable(Uniform(0, 1), :y)

        d = Model(df -> sqrt.(df.x .^ 2 + df.y .^ 2), :d)

        pf, _ = probability_of_failure(
            [d],
            df -> df.d .- 1,
            [x, y],
            MonteCarlo(1000000),
        )

        Random.seed!()

        @test round(4.0 * pf, digits = 3) == 3.141
    end

end
