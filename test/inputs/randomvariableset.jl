using Statistics, Random

rv1 = RandomVariable(Exponential(1), :x)
rv2 = RandomVariable(Exponential(1/2), :y)

rvs_ = [rv1 rv2]
corr_ = [1 0.8; 0.8 1]

@testset "RandomVariableSet" begin

    @testset "Constructor" begin
        @test isa(RandomVariableSet(rvs_), RandomVariableSet)
        @test isa(RandomVariableSet(rvs_, corr_), RandomVariableSet)
        @test isa(RandomVariableSet(members = rvs_), RandomVariableSet)
        @test isa(RandomVariableSet(members = rvs_,
            corr = corr_,
        ), RandomVariableSet)

        rvset = RandomVariableSet(rvs_, corr_)
        @test rvset.members == rvs_
        @test rvset.corr == corr_

        rvset = RandomVariableSet(members = rvs_)
        @test rvset.corr == [1 0; 0 1]

        @test_throws ErrorException("wrong dimension of correlation matrix") RandomVariableSet(rvs_,
            [1 0 0; 0 1 0; 0 0 1],
        )
    end

    @testset "sample" begin
        rvset = RandomVariableSet(rvs_, corr_)
        @test size(sample(rvset, 10)) == (10, 2)
        @test size(sample(rvset)) == (1, 2)
    end

    @testset "names" begin
        rvset = RandomVariableSet(rvs_, corr_)
        @test names(rvset) == [:x, :y]
    end

    @testset "mean" begin
        rvset = RandomVariableSet(rvs_, corr_)
        @test mean(rvset) == DataFrame(:x => 1.0, :y => 0.5)
    end

    @testset "to_standard_normal_space" begin
        rvset = RandomVariableSet(rvs_, corr_)

        Random.seed!(8128)

        samples = sample(rvset, 10^5)

        @test isapprox(Statistics.mean(samples.x), 1.0, atol = 0.01)
        @test isapprox(Statistics.mean(samples.y), 0.5, atol = 0.01)

        @test round(cor(samples.x, samples.y), digits = 2) == 0.77

        to_standard_normal_space!(rvset, samples)

        @test isapprox(abs(Statistics.mean(samples.x)), 0.0, atol = 0.01)
        @test isapprox(abs(Statistics.mean(samples.y)), 0.0, atol = 0.01)

        @test isapprox(std(samples.x), 1.0, atol = 0.01)
        @test isapprox(std(samples.y), 1.0, atol = 0.01)

        @test round(abs(cor(samples.x, samples.y)), digits = 2) == 0.0

        Random.seed!()
    end

end
