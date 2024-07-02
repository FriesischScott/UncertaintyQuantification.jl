@testset "Imprecise Probability of Failure" begin
    @testset "Double Loop" begin
        @testset "P-boxes only" begin
            X = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :X)
            Y = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :Y)

            pf = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], DoubleLoop(FORM())
            )

            @test pf.lb ≈ 0.0 atol = 1e-6
            @test pf.ub ≈ 0.006664164390407847 atol = 1e-6
        end
    end
end
