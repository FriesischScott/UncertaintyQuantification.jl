@testset "Imprecise Probability of Failure" begin
    @testset "Double Loop" begin
        @testset "P-boxes only" begin
            X = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :X)
            Y = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :Y)

            pf = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], DoubleLoop(FORM())
            )

            @test pf.lb ≈ cdf(Normal(2, sqrt(2)), -9) atol = 1e-6
            @test pf.ub ≈ cdf(Normal(-2, sqrt(8)), -9) atol = 1e-6
        end

        @testset "Interval - Distribution" begin

            X = Interval(-1, 1, :X)
            Y = RandomVariable(Normal(0, 2), :Y)

            pf = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], DoubleLoop(FORM())
            )

            @test pf.lb ≈ cdf(Normal(1, 2), -9) atol = 1e-6
            @test pf.ub ≈ cdf(Normal(-1, 2), -9) atol = 1e-6
        end

        @testset "Interval - p-box" begin
            
            X = Interval(-1, 1, :X)
            Y = ProbabilityBox{Normal}([Parameter(0, :μ), Interval(1, 2, :σ)], :Y)

            pf = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], DoubleLoop(FORM())
            )

            @test pf.lb ≈ cdf(Normal(1, 1), -9) atol = 1e-6
            @test pf.ub ≈ cdf(Normal(-1, 2), -9) atol = 1e-6
        end

    end
end
