@testset "Interval" begin
    name = :l
    lb = 0.14
    ub = 0.16
    @test_throws ErrorException(
        "lower bound parameter must be smaller than upper bound parameter for $name"
    ) Interval(ub, lb, name)
    interval = Interval(lb, ub, name)
    @test interval.lb == lb
    @test interval.ub == ub
    @test interval.name == name

    par = 0.13
    @test_throws ErrorException(
        "Choosen value $par is lower than Interval's lower bound $lb"
    ) UncertaintyQuantification.map_to_precise(par, interval)

    par = 0.17
    @test_throws ErrorException(
        "Choosen value $par is higher than Interval's upper bound $ub"
    ) UncertaintyQuantification.map_to_precise(par, interval)

    @test UncertaintyQuantification.map_to_precise(0.15, interval) ==
        Parameter(0.15, interval.name)

    interval = Interval(lb, ub, name)
    @test UncertaintyQuantification.sample(interval) == [lb, ub]
end