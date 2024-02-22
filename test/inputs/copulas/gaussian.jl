using Random

correlation = [1 0.7071; 0.7071 1]

@testset "GaussianCopula" begin
    copula = GaussianCopula(correlation)
    @test copula.correlation == correlation
    @test dimensions(copula) == 2

    @testset "sample" begin
        u = sample(copula, 10^5)

        @test corkendall(u[:, 1], u[:, 2]) ≈ 0.5 atol = 0.01
    end

    @testset "to_standard_normal_space" begin
        u = sample(copula, 10^5)
        s = to_standard_normal_space(copula, u)

        @test mean(s[:, 1]) ≈ 0 atol = 0.05
        @test median(s[:, 1]) ≈ 0 atol = 0.05
        @test std(s[:, 1]) ≈ 1 atol = 0.05

        @test mean(s[:, 2]) ≈ 0 atol = 0.05
        @test median(s[:, 2]) ≈ 0 atol = 0.05
        @test std(s[:, 2]) ≈ 1 atol = 0.05

        @test corkendall(s[:, 1], s[:, 2]) ≈ 0.0 atol = 0.01
    end

    @testset "to_copula_space" begin
        s = randn(10^5, 2)
        u = to_copula_space(copula, s)

        @test corkendall(u[:, 1], u[:, 2]) ≈ 0.5 atol = 0.01
    end
end
