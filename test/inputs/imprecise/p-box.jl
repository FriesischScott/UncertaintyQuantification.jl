@testset "P-box" begin
    name = :l

    params = [Interval(0.14, 0.16, :a), Interval(0.21, 0.23, :b), Interval(0, 1, :c)]
    @test_throws ErrorException(
        "Parameter mismatch for ProbabilityBox l: [:a, :b, :c] != [:a, :b]."
    ) ProbabilityBox{Uniform}(params, name)

    params = [Interval(0.14, 0.16, :a), Interval(0.21, 0.23, :b)]

    p_box = ProbabilityBox{Uniform}(params, name)
    @test p_box.parameters == params
    @test p_box.name == name

    @test UncertaintyQuantification.bounds(p_box) == ([0.14, 0.21], [0.16, 0.23])

    par = [0.13, 0.20]
    @test_throws ErrorException(
        "Values outside of parameter intervals for ProbabilityBox l."
    ) UncertaintyQuantification.map_to_precise(par, p_box)

    par = [0.17, 0.25]
    @test_throws ErrorException(
        "Values outside of parameter intervals for ProbabilityBox l."
    ) UncertaintyQuantification.map_to_precise(par, p_box)

    par = [0.15, 0.22]
    @test UncertaintyQuantification.map_to_precise(par, p_box) ==
        RandomVariable(Uniform(par...), p_box.name)

    p_box = ProbabilityBox{Normal}([Interval(0, 1, :μ), Interval(0.1, 1, :σ)], name)
    a, b = UncertaintyQuantification.sample(p_box, 0.25)
    @test a == quantile(Normal(0, 1), 0.25)
    @test b == quantile(Normal(1, 0.1), 0.25)

    a, b = UncertaintyQuantification.sample(p_box, 0.5)
    @test a == quantile(Normal(0, 0.1), 0.5)
    @test a == quantile(Normal(0, 1), 0.5)
    @test b == quantile(Normal(1, 0.1), 0.5)
    @test b == quantile(Normal(1, 1), 0.5)

    a, b = UncertaintyQuantification.sample(p_box, 0.75)
    @test a == quantile(Normal(0, 0.1), 0.75)
    @test b == quantile(Normal(1, 1), 0.75)

    p_box = ProbabilityBox{Normal}([Parameter(0, :μ), Interval(0.1, 1, :σ)], name)
    @test UncertaintyQuantification.bounds(p_box) == ([0.1], [1])

    a, b = UncertaintyQuantification.sample(p_box, 0.25)
    @test a == quantile(Normal(0, 1), 0.25)
    @test b == quantile(Normal(0, 0.1), 0.25)

    a, b = UncertaintyQuantification.sample(p_box, 0.5)
    @test a == 0.0
    @test b == 0.0

    a, b = UncertaintyQuantification.sample(p_box, 0.75)
    @test a == quantile(Normal(0, 0.1), 0.75)
    @test b == quantile(Normal(0, 1), 0.75)
end
