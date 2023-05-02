@testset "RandomVariable" begin
    dist = Normal(0, 1)
    name = :x
    x = RandomVariable(dist, name)

    @test isa(x, RandomVariable)
    @test x.dist == dist
    @test x.name == name
    @test dimensions(x) == 1

    @test size(sample(x)) == (1, 1)
    @test size(sample(x, 10)) == (10, 1)

    @test logpdf(x, -1.0) == logpdf(x.dist, -1.0)
    @test logpdf(x, 1.0) == logpdf(x.dist, 1.0)
    @test logpdf(x, 0.0) == logpdf(x.dist, 0.0)

    @test pdf(x, -1.0) == pdf(x.dist, -1.0)
    @test pdf(x, 1.0) == pdf(x.dist, 1.0)
    @test pdf(x, 0.0) == pdf(x.dist, 0.0)

    @test cdf(x, -1.0) == cdf(x.dist, -1.0)
    @test cdf(x, 1.0) == cdf(x.dist, 1.0)
    @test cdf(x, 0.0) == cdf(x.dist, 0.0)

    @test quantile(x, 0.0) == quantile(x.dist, 0.0)
    @test quantile(x, 0.5) == quantile(x.dist, 0.5)
    @test quantile(x, 1.0) == quantile(x.dist, 1.0)

    @test minimum(x) == minimum(x.dist)
    @test maximum(x) == maximum(x.dist)
    @test mean(x) == mean(x.dist)
    @test var(x) == var(x.dist)

    y = RandomVariable(Uniform(), :y)
    @test insupport(x, -1) == insupport(x.dist, -1)
    @test insupport(x, 0) == insupport(x.dist, 0)
    @test insupport(x, 1) == insupport(x.dist, 1)

    @testset "SNS Mappings" begin
        rvs = [
            RandomVariable(Normal(), :x),
            RandomVariable(Normal(3, 0.5), :x),
            RandomVariable(LogNormal(), :x),
            RandomVariable(Uniform(-2, 3), :x),
            RandomVariable(Exponential(0.5), :x),
        ]

        n = 10^5

        for rv in rvs
            samples = sample(rv, n)
            mapped = copy(samples)

            to_standard_normal_space!(rv, mapped)

            @test mean(mapped.x) ≈ 0.0 atol = 0.01
            @test std(mapped.x) ≈ 1.0 atol = 0.01

            to_physical_space!(rv, mapped)

            samples.x == mapped.x
        end
    end
end
