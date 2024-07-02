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
    @test UncertaintyQuantification.sample(interval) == [lb, ub]
end
