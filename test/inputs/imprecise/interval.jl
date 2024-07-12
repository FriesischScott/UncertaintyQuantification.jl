@testset "Interval" begin
    name = :l
    lb = 0.14
    ub = 0.16
    @test_throws ErrorException(
        "Lower bound parameter must be smaller than upper bound parameter for Interval $name.",
    ) Interval(ub, lb, name)
    interval = Interval(lb, ub, name)
    @test interval.lb == lb
    @test interval.ub == ub
    @test interval.name == name

    par = 0.13
    @test_throws ErrorException("0.13 not in [0.14, 0.16] for Interval l.") UncertaintyQuantification.map_to_precise(
        par, interval
    )

    par = 0.17
    @test_throws ErrorException("0.17 not in [0.14, 0.16] for Interval l.") UncertaintyQuantification.map_to_precise(
        par, interval
    )

    @test UncertaintyQuantification.map_to_precise(0.15, interval) ==
        Parameter(0.15, interval.name)

    interval = Interval(lb, ub, name)

    @test UncertaintyQuantification.sample(interval) == DataFrame(; l=interval)

    @testset "Interval propagation" begin
        X1 = ProbabilityBox{Normal}([Interval(-1, 2, :μ), Parameter(1, :σ)], :X1)
        X2 = ProbabilityBox{Normal}([Interval(-2, 1, :μ), Parameter(2, :σ)], :X2)
        X3 = RandomVariable(Normal(0, 1), :X3)
        X4 = Parameter(5, :X4)

        inputs = [X1, X2, X3, X4]
        models = Model(df -> df.X1 .^ 2 .+ df.X2 .+ df.X3 .+ df.X4, :g)

        df = sample(inputs, 500)
        evaluate!(models, df)

        @test eltype(df.g) == Interval
    end
end
