@testset "sobolindices" begin
    x = RandomVariable.(Uniform(-π, π), [:x1, :x2, :x3])

    ishigami1 = Model(
        df -> sin.(df.x1) + 7 .* sin.(df.x2) .^ 2 + 0.1 .* df.x3 .^ 4 .* sin.(df.x1), :f1
    )
    ishigami2 = Model(
        df -> sin.(df.x1) + 7 .* sin.(df.x2) .^ 2 + 0.05 .* df.x3 .^ 4 .* sin.(df.x1), :f2
    )

    n = 2200

    firstorder_analytical1 = [0.3138, 0.4424, 0.00]
    totaleffect_analytical1 = [0.5574, 0.4424, 0.2436]

    firstorder_analytical2 = [0.219, 0.687, 0.00]
    totaleffect_analytical2 = [0.3136, 0.687, 0.0946]

    @testset "Monte Carlo" begin
        Random.seed!(8128)

        si = sobolindices([ishigami1, ishigami2], x, [:f1, :f2], MonteCarlo(n))

        Random.seed!()

        @test all(isapprox(si[:f1].FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si[:f1].TotalEffect, totaleffect_analytical1; rtol=0.1))
        @test all(isapprox(si[:f2].FirstOrder, firstorder_analytical2; rtol=0.1))
        @test all(isapprox(si[:f2].TotalEffect, totaleffect_analytical2; rtol=0.1))
    end

    @testset "Sobol" begin
        Random.seed!(8128)

        si = sobolindices([ishigami1, ishigami2], x, [:f1, :f2], SobolSampling(n))

        @test all(isapprox(si[:f1].FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si[:f1].TotalEffect, totaleffect_analytical1; rtol=0.1))
    end

    @testset "Halton" begin
        Random.seed!(8128)

        si = sobolindices([ishigami1, ishigami2], x, [:f1, :f2], HaltonSampling(n))

        @test all(isapprox(si[:f1].FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si[:f1].TotalEffect, totaleffect_analytical1; rtol=0.1))
    end

    @testset "Latin Hypercube" begin
        Random.seed!(8128)

        si = sobolindices([ishigami1, ishigami2], x, [:f1, :f2], LatinHypercubeSampling(n))

        @test all(isapprox(si[:f1].FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si[:f1].TotalEffect, totaleffect_analytical1; rtol=0.1))
    end

    @testset "Polynomial Chaos Expansion" begin
        p = 6
        Ψ = PolynomialChaosBasis(LegendreBasis.([p, p, p]), p)

        gq = GaussQuadrature()
        pce, _ = polynomialchaos(x, ishigami1, Ψ, :f1, gq)

        si = sobolindices(pce)

        @test all(isapprox(si.FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical1; rtol=0.1))
    end
end
