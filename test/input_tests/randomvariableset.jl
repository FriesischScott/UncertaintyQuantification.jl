
using UncertaintyQuantification

using Distributions
using Test

rv1 = Normal(0,1);
rv2 = Normal(0,1);


rvs_ = [rv1 rv2]
names_ = ["Rv1" "Rv2"];
corr_ = [1 0.8; 0.8 1];

@testset "constructor tests" begin

    @test typeof(RandomVariableSet(rvs_, names_)) == RandomVariableSet;
    @test typeof(RandomVariableSet(rvs_, names_, corr_)) == RandomVariableSet;
    @test typeof(RandomVariableSet(members = rvs_, names = names_)) == RandomVariableSet;
    @test typeof(RandomVariableSet(members = rvs_, names = names_, corr = corr_)) == RandomVariableSet;

end