@testset "Imprecise Probability of Failure" begin
    @testset "Double Loop" begin
        @testset "P-boxes only" begin
            X = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :X)
            Y = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :Y)

            pf, x_lb, x_ub = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], DoubleLoop(FORM())
            )

            @test pf.lb ≈ cdf(Normal(2, sqrt(2)), -9) atol = 1e-6
            @test pf.ub ≈ cdf(Normal(-2, sqrt(8)), -9) atol = 1e-6
        end

        @testset "Interval - Distribution" begin
            X = Interval(-1, 1, :X)
            Y = RandomVariable(Normal(0, 2), :Y)

            pf, x_lb, x_ub = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], DoubleLoop(FORM())
            )

            @test pf.lb ≈ cdf(Normal(1, 2), -9) atol = 1e-6
            @test pf.ub ≈ cdf(Normal(-1, 2), -9) atol = 1e-6
        end

        @testset "P-box - Distribution" begin
            X = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :X)
            Y = RandomVariable(Normal(0, 2), :Y)

            pf, x_lb, x_ub = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], DoubleLoop(FORM())
            )

            @test pf.lb ≈ cdf(Normal(1, sqrt(5)), -9) atol = 1e-6
            @test pf.ub ≈ cdf(Normal(-1, sqrt(8)), -9) atol = 1e-6
        end

        @testset "Interval - p-box" begin
            X = Interval(-1, 1, :X)
            Y = ProbabilityBox{Normal}([Parameter(0, :μ), Interval(1, 2, :σ)], :Y)

            pf, x_lb, x_ub = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], DoubleLoop(FORM())
            )

            @test pf.lb ≈ cdf(Normal(1, 1), -9) atol = 1e-6
            @test pf.ub ≈ cdf(Normal(-1, 2), -9) atol = 1e-6

            @test x_lb ≈ [1, 1]
            @test x_ub ≈ [-1, 2]
        end
    end

    @testset "Random Slicing" begin
        @testset "P-boxes only" begin
            X = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :X)
            Y = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :Y)

            pf, _ = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], RandomSlicing(FORM())
            )

            pbox_analyt = ProbabilityBox{Normal}(
                [Interval(-2, 2, :μ), Interval(sqrt(2), sqrt(8), :σ)], :X
            )
            failure_analty = cdf(pbox_analyt, -9)

            @test pf.lb ≈ failure_analty.lb atol = 1e-6
            @test pf.ub ≈ failure_analty.ub atol = 1e-6
        end

        @testset "Interval - Distribution" begin
            X = Interval(-1, 1, :X)
            Y = RandomVariable(Normal(0, 2), :Y)

            pf, _ = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], RandomSlicing(FORM())
            )

            pbox_analyt = ProbabilityBox{Normal}(
                [Interval(-1, 1, :μ), Parameter(2, :σ)], :X
            )
            failure_analty = cdf(pbox_analyt, -9)

            @test pf.lb ≈ failure_analty.lb atol = 1e-6
            @test pf.ub ≈ failure_analty.ub atol = 1e-6
        end

        @testset "P-box - Distribution" begin
            X = ProbabilityBox{Normal}([Interval(-1, 1, :μ), Interval(1, 2, :σ)], :X)
            Y = RandomVariable(Normal(0, 2), :Y)

            pf, _ = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], RandomSlicing(FORM())
            )

            pbox_analyt = ProbabilityBox{Normal}(
                [Interval(-1, 1, :μ), Interval(sqrt(5), sqrt(8), :σ)], :X
            )
            failure_analty = cdf(pbox_analyt, -9)

            @test pf.lb ≈ cdf(Normal(1, sqrt(5)), -9) atol = 1e-6
            @test pf.ub ≈ cdf(Normal(-1, sqrt(8)), -9) atol = 1e-6
        end

        @testset "Interval - p-box" begin
            X = Interval(-1, 1, :X)
            Y = ProbabilityBox{Normal}([Parameter(0, :μ), Interval(1, 2, :σ)], :Y)

            pf, _ = probability_of_failure(
                UQModel[], df -> 9 .+ df.X .+ df.Y, [X, Y], RandomSlicing(FORM())
            )

            pbox_analyt = ProbabilityBox{Normal}(
                [Interval(-1, 1, :μ), Interval(1, 2, :σ)], :X
            )
            failure_analty = cdf(pbox_analyt, -9)

            @test pf.lb ≈ failure_analty.lb atol = 1e-6
            @test pf.ub ≈ failure_analty.ub atol = 1e-6
        end
    end
end
