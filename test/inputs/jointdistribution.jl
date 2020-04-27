using Random, DataFrames

rv1 = RandomVariable(Exponential(1), :x)
rv2 = RandomVariable(Exponential(1/2), :y)

marginals = [rv1 rv2]
copula = GaussianCopula([1 0.8; 0.8 1])

@testset "JointDistribution" begin

    @testset "Constructor" begin
        @test isa(JointDistribution(marginals, copula), JointDistribution)

        jd = JointDistribution(marginals, copula)
        @test jd.marginals == marginals
        @test jd.copula == copula

        @test_throws ErrorException("Dimension mismatch between copula and marginals") JointDistribution(marginals,
            GaussianCopula([1 0 0; 0 1 0; 0 0 1]),
        )
    end

    @testset "sample" begin
        jd = JointDistribution(marginals, copula)
        @test size(sample(jd, 10)) == (10, 2)
        @test size(sample(jd)) == (1, 2)
    end

    @testset "names" begin
        jd = JointDistribution(marginals, copula)
        @test names(jd) == [:x, :y]
    end

    @testset "mean" begin
        jd = JointDistribution(marginals, copula)
        @test mean(jd) == DataFrame(:x => 1.0, :y => 0.5)
    end

    @testset "to_standard_normal_space" begin
        jd = JointDistribution(marginals, copula)

        Random.seed!(8128)

        samples = sample(jd, 10^5)

        @test isapprox(mean(samples.x), 1.0, atol = 0.01)
        @test isapprox(mean(samples.y), 0.5, atol = 0.01)

        @test round(cor(samples.x, samples.y), digits = 2) == 0.77

        to_standard_normal_space!(jd, samples)

        @test isapprox(abs(mean(samples.x)), 0.0, atol = 0.01)
        @test isapprox(abs(mean(samples.y)), 0.0, atol = 0.01)

        @test isapprox(std(samples.x), 1.0, atol = 0.01)
        @test isapprox(std(samples.y), 1.0, atol = 0.01)

        @test round(abs(cor(samples.x, samples.y)), digits = 2) == 0.0

        Random.seed!()
    end

    @testset "to_physical_space" begin
    jd = JointDistribution(marginals, copula)

    Random.seed!(8128)

    samples = DataFrame(:x => rand(Normal(), 10^5), :y => rand(Normal(), 10^5))

    to_physical_space!(jd, samples)

    @test isapprox(mean(samples.x), 1.0, atol = 0.01)
    @test isapprox(mean(samples.y), 0.5, atol = 0.01)

    @test round(cor(samples.x, samples.y), digits = 2) == 0.77

    Random.seed!()
end

end
