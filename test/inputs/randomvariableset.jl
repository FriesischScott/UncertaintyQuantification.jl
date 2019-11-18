rv1 = RandomVariable(Normal(0, 1), :Rv1)
rv2 = RandomVariable(Normal(0, 1), :RV2)

rvs_ = [rv1 rv2]
corr_ = [1 0.8; 0.8 1]

@testset "RandomVariableSet" begin

    @testset "Constructor" begin
        @test isa(RandomVariableSet(rvs_), RandomVariableSet)
        @test isa(RandomVariableSet(rvs_, corr_), RandomVariableSet)
        @test isa(RandomVariableSet(members = rvs_), RandomVariableSet)
        @test isa(RandomVariableSet(
            members = rvs_,
            corr = corr_,
        ), RandomVariableSet)

        rvset = RandomVariableSet(rvs_, corr_)
        @test rvset.members == rvs_
        @test rvset.corr == corr_

        rvset = RandomVariableSet(members = rvs_)
        @test rvset.corr == [1 0; 0 1]

        @test_throws ErrorException("wrong dimension of correlation matrix") RandomVariableSet(
            rvs_,
            [1 0 0; 0 1 0; 0 0 1],
        )

    end

    @testset "sample" begin
        rvset = RandomVariableSet(rvs_, corr_)
        @test size(sample(rvset, 10)) == (10, 2)
        @test size(sample(rvset)) == (1, 2)
    end
end
