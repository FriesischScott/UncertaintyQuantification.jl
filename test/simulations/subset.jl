@testset "SubSetSimulation" begin
    proposal = Normal()
    subset = SubSetSimulation(2000, 0.2, 10, proposal)

    @test isa(subset, SubSetSimulation)
    @test subset.n == 2000
    @test subset.target == 0.2
    @test subset.levels == 10
    @test subset.proposal == proposal

    @test_throws ErrorException("proposal must be a symmetric distribution") SubSetSimulation(
        2000, 0.2, 10, Exponential()
    )
    @test_throws ErrorException("proposal must be centered in 0") SubSetSimulation(
        2000, 0.2, 10, Uniform()
    )
    @test_logs (:warn, "A proposal pdf with large variance (≥ 2) can be inefficient.") SubSetSimulation(
        2000, 0.2, 10, Uniform(-4, 4)
    )
end

@testset "SubSetInfinity" begin
    subset = SubSetInfinity(2000, 0.2, 10, 0.5)

    @test isa(subset, SubSetInfinity)
    @test subset.n == 2000
    @test subset.target == 0.2
    @test subset.levels == 10
    @test subset.s == 0.5

    @test_throws ErrorException("standard deviation must be between 0.0 and 1.0") SubSetInfinity(
        2000, 0.2, 10, 2.0
    )
    @test_throws ErrorException("standard deviation must be between 0.0 and 1.0") SubSetInfinity(
        2000, 0.2, 10, -1.0
    )
end

@testset "SubSetInfinityAdaptive" begin
    subset = SubSetInfinityAdaptive(2000, 0.2, 10, 4, 1, 0.5)

    @test isa(subset, SubSetInfinityAdaptive)
    @test subset.n == 2000
    @test subset.target == 0.2
    @test subset.levels == 10
    @test subset.Na == 4
    @test subset.λ == 1
    @test subset.s == 0.5

    subset = SubSetInfinityAdaptive(2000, 0.2, 10, 4, 1)

    @test isa(subset, SubSetInfinityAdaptive)
    @test subset.n == 2000
    @test subset.target == 0.2
    @test subset.levels == 10
    @test subset.Na == 4
    @test subset.λ == 1
    @test subset.s == 1

    subset = SubSetInfinityAdaptive(2000, 0.2, 10, 4)

    @test isa(subset, SubSetInfinityAdaptive)
    @test subset.n == 2000
    @test subset.target == 0.2
    @test subset.levels == 10
    @test subset.Na == 4
    @test subset.λ == 1
    @test subset.s == 1


    @test_throws ErrorException(
        "standard deviation must be between 0.0 and 1.0"
    ) SubSetInfinityAdaptive(2000, 0.1, 10, 2, 1, 3)

    @test_throws ErrorException(
        "Scaling parameter must be between 0.0 and 1.0. A good initial choice is 1.0"
    ) SubSetInfinityAdaptive(2000, 0.1, 10, 4, 2.0)

    @test_throws ErrorException(
        "Scaling parameter must be between 0.0 and 1.0. A good initial choice is 1.0"
    ) SubSetInfinityAdaptive(2000, 0.1, 10, 4, -2.0)
    
    @test_throws ErrorException("Number of partitions Na must be a multiple of `n` * `target`") SubSetInfinityAdaptive(
        2000, 0.1, 10, 9, 1
    )

    @test_throws ErrorException("Number of partitions Na must be less than `n` * `target`") SubSetInfinityAdaptive(
        2000, 0.1, 10, 400, 1
    )
end
