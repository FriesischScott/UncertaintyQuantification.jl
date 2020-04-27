using Random

correlation = [1 0.8; 0.8 1]

@testset "GaussianCopula" begin

    copula = GaussianCopula(correlation)
    @test copula.correlation == correlation
    @test dimension(copula) == 2
    
    @testset "sample" begin
        Random.seed!(8128)

        u = sample(copula, 10^5)

        @test round(cor(u[:, 1], u[:, 2]), digits = 2) == 0.79

        Random.seed!()
    end
    
    @testset "to_standard_normal_space" begin
        Random.seed!(8128)

        u = sample(copula, 10^5)
        s = to_standard_normal_space(copula, u)

        @test isapprox(abs(mean(s[:, 1])), 0.0, atol = 0.01)
        @test isapprox(abs(mean(s[:, 2])), 0.0, atol = 0.01)

        @test isapprox(std(s[:, 1]), 1.0, atol = 0.01)
        @test isapprox(std(s[:, 2]), 1.0, atol = 0.01)

        @test round(abs(cor(s[:, 1], s[:, 2])), digits = 2) == 0.0

        Random.seed!()
    end

    @testset "to_copula_space" begin
        Random.seed!(8128)

        s = randn(10^5, 2)
        u = to_copula_space(copula, s)

        @test round(cor(u[:, 1], u[:, 2]), digits = 2) == 0.79


        Random.seed!()
    end
end