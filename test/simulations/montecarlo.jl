@testset "MonteCarlo" begin
    mc = MonteCarlo(1000)

    @test isa(mc, MonteCarlo)
    @test mc.n == 1000
end

@testset "QuasiMonteCarlo" begin

    # Make sure that all QMC Types implement qmc_samples
    for qmc in subtypes(AbstractQuasiMonteCarlo)
        @test hasmethod(qmc_samples, (qmc, Integer))
    end

    @testset "SobolSampling" begin
        sobol = SobolSampling(1000)

        @test isa(sobol, SobolSampling)
        @test sobol.n == 1000

        @testset "sample" begin
            inputs = [RandomVariable.(Uniform(), [:a, :b]); Parameter(1, :c)]

            samples = sample(inputs, SobolSampling(4))

            @test isapprox(samples.a, [0.375, 0.875, 0.625, 0.125])
            @test isapprox(samples.b, [0.375, 0.875, 0.125, 0.625])
            @test samples.c == [1.0, 1.0, 1.0, 1.0]
        end
    end

    @testset "HaltonSampling" begin
        halton = HaltonSampling(1000)

        @test isa(halton, HaltonSampling)
        @test halton.n == 1000

        @testset "sample" begin
            inputs = [
                RandomVariable.([Uniform(-1, 0), Uniform()], [:a, :b])
                Parameter(1, :c)
            ]

            samples = sample(inputs, halton)

            @test mean(samples.a) ≈ -0.5 rtol = 0.05
            @test mean(samples.b) ≈ 0.5 rtol = 0.05
            @test samples.c == ones(halton.n)
        end
    end

    @testset "LatticeRuleSampling" begin
        lattice = LatticeRuleSampling(1000)

        @test isa(lattice, LatticeRuleSampling)
        @test lattice.n == 1000

        @testset "sample" begin
            inputs = [RandomVariable.(Uniform(), [:a, :b]); Parameter(1, :c)]

            samples = sample(inputs, lattice)
            @test mean(samples.a) ≈ 0.5 rtol = 0.05
            @test mean(samples.b) ≈ 0.5 rtol = 0.05
            @test samples.c == ones(lattice.n)
        end
    end

    @testset "LatinHypercubeSampling" begin
        lhs = LatinHypercubeSampling(1000)

        inputs = RandomVariable.([Uniform(-1, 1), Uniform(0, 3)], [:a, :b])

        samples = sample(inputs, lhs)

        h = fit(Histogram, samples.a, range(-1, 1, 1001))
        @test all(h.weights .== 1)

        h = fit(Histogram, samples.b, range(0, 3, 1001))
        @test all(h.weights .== 1)
    end
end
