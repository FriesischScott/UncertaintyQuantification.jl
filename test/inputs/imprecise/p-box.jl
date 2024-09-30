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

    @testset "Quantile and sampling" begin
        name = :l

        p_box = ProbabilityBox{Normal}([Interval(0, 1, :μ), Interval(0.1, 1, :σ)], name)
        a = UncertaintyQuantification.quantile(p_box, 0.25)
        @test a.lb == quantile(Normal(0, 1), 0.25)
        @test a.ub == quantile(Normal(1, 0.1), 0.25)

        a = UncertaintyQuantification.quantile(p_box, 0.5)
        @test a.lb == quantile(Normal(0, 0.1), 0.5)
        @test a.lb == quantile(Normal(0, 1), 0.5)
        @test a.ub == quantile(Normal(1, 0.1), 0.5)
        @test a.ub == quantile(Normal(1, 1), 0.5)

        a = UncertaintyQuantification.quantile(p_box, 0.75)
        @test a.lb == quantile(Normal(0, 0.1), 0.75)
        @test a.ub == quantile(Normal(1, 1), 0.75)

        p_box = ProbabilityBox{Normal}([Parameter(0, :μ), Interval(0.1, 1, :σ)], name)
        @test UncertaintyQuantification.bounds(p_box) == ([0.1], [1])

        a = UncertaintyQuantification.quantile(p_box, 0.25)
        @test a.lb == quantile(Normal(0, 1), 0.25)
        @test a.ub == quantile(Normal(0, 0.1), 0.25)

        a = UncertaintyQuantification.quantile(p_box, 0.5)
        @test a == Interval(0.0, 0.0, :l)

        a = UncertaintyQuantification.quantile(p_box, 0.75)
        @test a.lb == quantile(Normal(0, 0.1), 0.75)
        @test a.ub == quantile(Normal(0, 1), 0.75)
    end

    @testset "Inverse quantile and transformations" begin
        p_box = ProbabilityBox{Uniform}([Interval(0, 0.2, :a), Interval(0.5, 1, :b)], name)

        Nsamples = 1000

        u = rand(Nsamples)
        x = quantile.(Ref(p_box), u)

        u_back = UncertaintyQuantification.reverse_quantile.(Ref(p_box), x)

        @test all(abs.(u_back .- u) .<= 10^-10)

        SNS_distribution = RandomVariable(Normal(0, 1), name)
        SNS_samples = sample(SNS_distribution, Nsamples)

        SNS_samples_before = deepcopy(SNS_samples)

        to_physical_space!(p_box, SNS_samples)
        to_standard_normal_space!(p_box, SNS_samples)

        @test all(abs.(SNS_samples[!, :l] .- SNS_samples_before[!, :l]) .<= 10^-10)

        p_box = ProbabilityBox{Uniform}([Interval(0, 0.2, :a), Parameter(1, :b)], name)

        u = rand(Nsamples)
        x = quantile.(Ref(p_box), u)

        u_back = UncertaintyQuantification.reverse_quantile.(Ref(p_box), x)

        @test all(abs.(u_back .- u) .<= 10^-10)

        SNS_distribution = RandomVariable(Normal(0, 1), name)
        SNS_samples = sample(SNS_distribution, Nsamples)

        SNS_samples_before = deepcopy(SNS_samples)

        to_physical_space!(p_box, SNS_samples)
        to_standard_normal_space!(p_box, SNS_samples)

        @test all(abs.(SNS_samples[!, :l] .- SNS_samples_before[!, :l]) .<= 10^-10)
    end
end
