@testset "P-box" begin
    name = :l
    lb = [0.14, 0.21]
    ub = [0.16, 0.23]
    dist = x -> Uniform(x...)

    @test_throws ErrorException(
        "lower bound parameters must be smaller than upper bound parameters for $name"
    ) ProbabilityBox(ub, lb,dist, name)

    p_box = ProbabilityBox(lb, ub, dist, name)
    @test p_box.lb == lb
    @test p_box.ub == ub
    @test p_box.dist == dist
    @test p_box.name == name

    par = [0.13, 0.20, 0,21]
    @test_throws ErrorException(
        "number of parameters $par must be equals to the number of parameter need by $name",
    ) UncertaintyQuantification.map_to_precise(par, p_box)

    par = [0.13, 0.20]
    @test_throws ErrorException(
        "One or more values in $par are lower than p-box's lower bound $lb",
    ) UncertaintyQuantification.map_to_precise(par, p_box)

    par = [0.17, 0.25]
    @test_throws ErrorException(
        "One or more values in $par are higher than p-box's upper bound $ub",
    ) UncertaintyQuantification.map_to_precise(par, p_box)

    par = [0.15, 0.22]
    @test UncertaintyQuantification.map_to_precise(par, p_box) ==
        RandomVariable(Uniform(par...), p_box.name)
end