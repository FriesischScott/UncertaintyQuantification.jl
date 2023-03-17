@testset "sobolindices" begin
    x = RandomVariable.(Uniform(-π, π), [:x1, :x2, :x3])

    ishigami1 = Model(
        df -> sin.(df.x1) + 7 .* sin.(df.x2) .^ 2 + 0.1 .* df.x3 .^ 4 .* sin.(df.x1), :f1
    )
    ishigami2 = Model(
        df -> sin.(df.x1) + 7 .* sin.(df.x2) .^ 2 + 0.05 .* df.x3 .^ 4 .* sin.(df.x1), :f2
    )

    n_mc = 4000
    n_qmc = 2300

    firstorder_analytical1 = [0.3138, 0.4424, 0.00]
    totaleffect_analytical1 = [0.5574, 0.4424, 0.2436]

    firstorder_analytical2 = [0.219, 0.687, 0.00]
    totaleffect_analytical2 = [0.3136, 0.687, 0.0946]

    @testset "Monte Carlo" begin
        Random.seed!(8128)

        si = sobolindices(ishigami1, x, :f1, MonteCarlo(n_mc))

        Random.seed!()

        @test all(isapprox(si.FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical1; rtol=0.1))
    end

    @testset "Sobol" begin
        si = sobolindices(ishigami1, x, :f1, SobolSampling(n_qmc))

        @test all(isapprox(si.FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical1; rtol=0.1))
    end

    @testset "Sobol Nultiple Outputs" begin
        si = sobolindices([ishigami1, ishigami2], x, [:f1, :f2], SobolSampling(n_qmc))

        @test all(isapprox(si[:f1].FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si[:f1].TotalEffect, totaleffect_analytical1; rtol=0.1))
        @test all(isapprox(si[:f2].FirstOrder, firstorder_analytical2; rtol=0.1))
        @test all(isapprox(si[:f2].TotalEffect, totaleffect_analytical2; rtol=0.1))
    end

    @testset "Halton" begin
        si = sobolindices(ishigami1, x, :f1, HaltonSampling(n_qmc))

        @test all(isapprox(si.FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical1; rtol=0.1))
    end

    @testset "Latin Hypercube" begin
        si = sobolindices(ishigami1, x, :f1, LatinHypercubeSampling(n_qmc))

        @test all(isapprox(si.FirstOrder, firstorder_analytical1; rtol=0.1))
        @test all(isapprox(si.TotalEffect, totaleffect_analytical1; rtol=0.1))
    end

    @testset "Convenience Functions" begin
        x1 = RandomVariable(Uniform(-π, +π), :x1)
        x2 = RandomVariable(Uniform(0, +π), :x2)
        t1 = Model(df -> sin.(df.x1), :t1)
        t2 = Model(df -> cos.(df.x1), :t2)

        @test hasmethod(sobolindices, typeof.([[t1; t2], x1, :t1, SobolSampling(n_mc)]))
        @test hasmethod(sobolindices, typeof.([t1, [x1; x2], :t1, SobolSampling(n_mc)]))
        @test hasmethod(sobolindices, typeof.([t1, x1, [:t1; :t2], SobolSampling(n_mc)]))
        @test hasmethod(
            sobolindices, typeof.([[t1; t2], [x1; x2], :t1, SobolSampling(n_mc)])
        )
        @test hasmethod(
            sobolindices, typeof.([[t1; t2], x1, [:t1, :t2], SobolSampling(n_mc)])
        )
        @test hasmethod(
            sobolindices, typeof.([t1, [x1; x2], [:t1; :t2], SobolSampling(n_mc)])
        )
        @test hasmethod(sobolindices, typeof.([t1, x1, :t1, SobolSampling(n_mc)]))
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