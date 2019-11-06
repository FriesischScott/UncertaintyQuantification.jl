@testset "Probability of Failure" begin

    @testset "Monte Carlo" begin
        Random.seed!(8128)

        a = RandomVariable(Uniform(1, 10), "a")
        b = RandomVariable(Normal(2, 1), "b")

        model = Model(x -> x.a .^ 2 .+ x.b .^ 3, "y")

        pf, _ = probabilityOfFailure(
            [model],
            x -> x.y,
            [a, b],
            MonteCarlo(100000)
        )

        Random.seed!()

        @test pf == 3.0e-5
    end

end
