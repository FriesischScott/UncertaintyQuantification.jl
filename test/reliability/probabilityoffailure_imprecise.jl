@testset "Imprecise Probability of Failure" begin
    x1 = ProbabilityBox{Normal}([Interval(0.5, 1.5, :μ), Parameter(1, :σ)], :x1)
    x2 = ProbabilityBox{Normal}([Interval(0.5, 1.5, :μ), Parameter(1, :σ)], :x2)

    g = @. Model(df -> df.x1 + df.x2, :g)

    F = df -> df.g
    @testset "Double Loop" begin
        pf = probability_of_failure(g, F, [x1, x2], MonteCarlo(1_000_000))

        @test pf.lb ≈ cdf(Normal(), -sqrt(2) * 1.5) atol = 0.01
        @test pf.ub ≈ cdf(Normal(), -sqrt(2) * 0.5) atol = 0.01
    end

    @testset "Interval Monte Carlo" begin
        pf = probability_of_failure(g, F, [x1, x2], MonteCarlo(10_000))

        @test pf.lb ≈ cdf(Normal(), -sqrt(2) * 1.5) atol = 0.01
        @test pf.ub ≈ cdf(Normal(), -sqrt(2) * 0.5) atol = 0.01
    end
end
