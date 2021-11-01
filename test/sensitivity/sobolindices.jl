@testset "sobolindices" begin
    x = RandomVariable.(Uniform(-π, π), [:x1, :x2, :x3])

    ishigami = Model(
        df -> sin.(df.x1) + 7 .* sin.(df.x2) .^ 2 + 0.1 .* df.x3 .^ 4 .* sin.(df.x1), :f
    )

    n = 1000

    firstorder_analytical = [0.3138, 0.4424, 0.00]
    totaleffect_analytical = [0.5574, 0.4424, 0.2436]

    @testset "Monte Carlo" begin
        Random.seed!(8128)

        si = sobolindices([ishigami], x, :f, MonteCarlo(n))

        Random.seed!()

        @test all(isapprox(si.FirstOrder, firstorder_analytical; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical; rtol=0.1))
    end

    @testset "Sobol" begin
        Random.seed!(8128)

        si = sobolindices([ishigami], x, :f, SobolSampling(n))

        Random.seed!()

        @test all(isapprox(si.FirstOrder, firstorder_analytical; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical; rtol=0.1))
    end

    @testset "Halton" begin
        Random.seed!(8128)

        si = sobolindices([ishigami], x, :f, MonteCarlo(n))

        Random.seed!()

        @test all(isapprox(si.FirstOrder, firstorder_analytical; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical; rtol=0.1))
    end

    @testset "Latin Hypercube" begin
        Random.seed!(8128)

        si = sobolindices([ishigami], x, :f, MonteCarlo(n))

        Random.seed!()

        @test all(isapprox(si.FirstOrder, firstorder_analytical; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical; rtol=0.1))
    end
end
