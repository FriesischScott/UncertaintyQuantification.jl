rv1 = RandomVariable(Normal(0, 1), "Rv1")
rv2 = RandomVariable(Normal(0, 1), "RV2")

rvs_ = [rv1 rv2]
corr_ = [1 0.8; 0.8 1]

@testset "RandomVariableSet" begin

    @test typeof(RandomVariableSet(rvs_)) == RandomVariableSet
    @test typeof(RandomVariableSet(rvs_, corr_)) == RandomVariableSet
    @test typeof(RandomVariableSet(members = rvs_)) == RandomVariableSet
    @test typeof(RandomVariableSet(
        members = rvs_,
        corr = corr_,
    )) == RandomVariableSet

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
