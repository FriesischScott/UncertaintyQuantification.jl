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

        @test abs(mean(s[:, 1])) ≈ 0.0 atol = 0.01
        @test abs(mean(s[:, 2])) ≈ 0.0 atol = 0.01

        @test std(s[:, 1]) ≈ 1.0 atol = 0.01
        @test std(s[:, 2]) ≈ 1.0 atol = 0.01

        @test cor(s[:, 1], s[:, 2]) ≈ 0.00 atol = 0.01
    end

    @testset "to_copula_space" begin
        s = randn(10^5, 2)
        u = to_copula_space(copula, s)

        @test cor(u[:, 1], u[:, 2]) ≈ 0.79 atol = 0.01
    end
end
