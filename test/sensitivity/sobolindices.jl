@testset "sobolindices" begin

    x = RandomVariable.(Uniform(-π, π), [:x1, :x2, :x3])

    ishigami = Model(
        df -> sin.(df.x1) + 7 .* sin.(df.x2) .^ 2 + 0.1 .* df.x3 .^ 4 .* sin.(df.x1),
        :f,
    )

    n = 1000

    @testset "Monte Carlo" begin
        Random.seed!(8128)

        si = sobolindices([ishigami], x, :f, MonteCarlo(n))

        Random.seed!()

        @test isapprox(si.FirstOrder, [0.314, 0.442, 0.00], rtol = 0.1) |> all
    end

    @testset "Sobol" begin
        Random.seed!(8128)

        si = sobolindices([ishigami], x, :f, SobolSampling(n))

        Random.seed!()

        @test isapprox(si.FirstOrder, [0.314, 0.442, 0.00], rtol = 0.1) |> all
    end

    @testset "Halton" begin
        Random.seed!(8128)

        si = sobolindices([ishigami], x, :f, MonteCarlo(n))

        Random.seed!()

        @test isapprox(si.FirstOrder, [0.314, 0.442, 0.00], rtol = 0.1) |> all
    end

    @testset "Latin Hypercube" begin
        Random.seed!(8128)

        si = sobolindices([ishigami], x, :f, MonteCarlo(n))

        Random.seed!()

        @test isapprox(si.FirstOrder, [0.314, 0.442, 0.00], rtol = 0.1) |> all
    end
end
