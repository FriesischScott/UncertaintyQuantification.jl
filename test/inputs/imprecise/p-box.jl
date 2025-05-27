@testset "P-box" begin
    params = Dict(
        :a => Interval(0.14, 0.16), :b => Interval(0.21, 0.23), :c => Interval(0, 1)
    )
    @test_throws ErrorException(
        "Parameter mismatch for ProbabilityBox [:a, :b, :c] != [:a, :b]."
    ) ProbabilityBox{Uniform}(params)

    params = Dict(:a => Interval(0.14, 0.16), :b => Interval(0.21, 0.23))

    p_box = ProbabilityBox{Uniform}(params)
    @test p_box.parameters == params
    @test p_box.lb == 0.14
    @test p_box.ub == 0.23

    @test UncertaintyQuantification.bounds(p_box) == ([0.14, 0.21], [0.16, 0.23])

    par = [0.13, 0.20]
    @test_throws ErrorException("Values outside of parameter intervals for ProbabilityBox") UncertaintyQuantification.map_to_distribution(
        par, p_box
    )

    par = [0.17, 0.25]
    @test_throws ErrorException("Values outside of parameter intervals for ProbabilityBox") UncertaintyQuantification.map_to_distribution(
        par, p_box
    )

    par = [0.15, 0.22]
    @test UncertaintyQuantification.map_to_distribution(par, p_box) == Uniform(par...)

    @testset "Quantile and sampling" begin
        p_box = ProbabilityBox{Normal}(Dict(:μ => Interval(0, 1), :σ => Interval(0.1, 1)))
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

        p_box = ProbabilityBox{Normal}(Dict(:μ => 0, :σ => Interval(0.1, 1)))
        @test UncertaintyQuantification.bounds(p_box) == ([0.1], [1])

        a = UncertaintyQuantification.quantile(p_box, 0.25)
        @test a.lb == quantile(Normal(0, 1), 0.25)
        @test a.ub == quantile(Normal(0, 0.1), 0.25)

        a = UncertaintyQuantification.quantile(p_box, 0.5)
        @test a == Interval(0.0, 0.0)

        a = UncertaintyQuantification.quantile(p_box, 0.75)
        @test a.lb == quantile(Normal(0, 0.1), 0.75)
        @test a.ub == quantile(Normal(0, 1), 0.75)
    end

    @testset "Inverse quantile" begin
        p_box = ProbabilityBox{Uniform}(
            Dict(:a => Interval(0, 0.2), :b => Interval(0.5, 1))
        )

        Nsamples = 1000

        u = rand(Nsamples)
        x = quantile.(Ref(p_box), u)

        u_back = UncertaintyQuantification.reverse_quantile.(Ref(p_box), x)

        @test all(abs.(u_back .- u) .<= 10^-10)
    end
end
