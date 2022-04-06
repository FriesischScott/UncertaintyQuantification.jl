@testset "PolynomialChaosExpansion" begin
    x = RandomVariable.(Uniform(-1, 1), [:x1, :x2])

    model = Model(df -> begin
        π .* (df.x1 .- 1) .* sin.(π .* df.x1) .* (1 .- df.x2 .^ 2)
    end, :y)

    p = 8
    Ψ = PolynomialChaosBasis(LegendreBasis.([p, p]), p)

    @testset "LeastSquares" begin
        ls = LeastSquares(SobolSampling(1000))
        pce, _, _ = polynomialchaos(x, model, Ψ, :y, ls)

        @test mean(pce) ≈ 2 / 3 rtol = 0.001
        @test var(pce) ≈ 2.931407384059896 rtol = 0.001
    end

    @testset "GaussQuadrature" begin
        gq = GaussQuadrature()
        pce, _ = polynomialchaos(x, model, Ψ, :y, gq)
        @test mean(pce) ≈ 2 / 3
        @test var(pce) ≈ 2.9313966486926053
    end

    @testset "evaluate" begin
        gq = GaussQuadrature()
        pce, samples = polynomialchaos(x, model, Ψ, :y, gq)

        data = copy(samples)
        evaluate!(pce, data)

        @test sum((samples.y .- data.y) .^ 2) ≈ 0.0 atol = 0.01
    end

    @testset "sample" begin
        gq = GaussQuadrature()
        pce, _ = polynomialchaos(x, model, Ψ, :y, gq)

        samples = sample(pce, 100)
        data = copy(samples)

        evaluate!(model, data)

        @test sum((samples.y .- data.y) .^ 2) ≈ 0.0 atol = 0.01
    end
end
