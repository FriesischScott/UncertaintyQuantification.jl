@testset "SubSetSimulation" begin
    proposal = Normal()
    subset = SubSetSimulation(2000, 0.2, 10, proposal)

    @test isa(subset, SubSetSimulation)
    @test subset.n == 2000
    @test subset.target == 0.2
    @test subset.levels == 10
    @test subset.proposal == proposal

    @test_throws ErrorException("proposal must be a symmetric distribution") SubSetSimulation(2000, 0.2, 10, Exponential())
end