@testset "PolynomialChaosExpansion" begin
    x = RandomVariable.([Uniform(-2, 0), Normal(-1, 0.5), Uniform(0, 1)], [:x1, :x2, :x3])

    model = Model(df -> begin
        return df.x1 .+ df.x2 .* df.x3
    end, :y)

    p = 8
    Ψ = PolynomialChaosBasis([LegendreBasis(p), HermiteBasis(), LegendreBasis(p)], p)

    @testset "LeastSquares" begin
        ls = LeastSquares(SobolSampling(1000))
        pce, _, _ = polynomialchaos(x, model, Ψ, :y, ls)

        @test mean(pce) ≈ -1.5 rtol = 1e-10
        @test var(pce) ≈ 0.5 rtol = 1e-10
    end

    @testset "GaussQuadrature" begin
        gq = GaussQuadrature()
        pce, _ = polynomialchaos(x, model, Ψ, :y, gq)

        @test mean(pce) ≈ -1.5 rtol = 1e-10
        @test var(pce) ≈ 0.5 rtol = 1e-10
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
