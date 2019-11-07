@testset "RandomVariable" begin

    dist = Normal(0, 1)
    name = "x"
    x = RandomVariable(dist, name)

    @test isa(x, RandomVariable)
    @test x.dist == dist
    @test x.name == name

    @test size(sample(x)) == (1, 1)
    @test size(sample(x, 10)) == (10, 1)
end
