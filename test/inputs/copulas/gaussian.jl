using Random

correlation = [1 0.8; 0.8 1]

@testset "GaussianCopula" begin
    copula = GaussianCopula(correlation)
    @test copula.correlation == correlation
    @test dimensions(copula) == 2

    @testset "sample" begin
        u = sample(copula, 10^5)

        @test cor(u[:, 1], u[:, 2]) ≈ 0.79 atol = 0.01
    end

    @testset "to_standard_normal_space" begin
        u = sample(copula, 10^5)
        s = to_standard_normal_space(copula, u)

        h0 = OneSampleADTest(s[:, 1], Normal())
        @test pvalue(h0) > 0.05

        h0 = OneSampleADTest(s[:, 2], Normal())
        @test pvalue(h0) > 0.05

        h0 = CorrelationTest(s[:, 1], s[:, 2])
        @test pvalue(h0) > 0.05
    end

    @testset "to_copula_space" begin
        s = randn(10^5, 2)
        u = to_copula_space(copula, s)

        @test cor(u[:, 1], u[:, 2]) ≈ 0.79 atol = 0.01
    end
end
