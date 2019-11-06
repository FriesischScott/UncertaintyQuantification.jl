@testset "RandomVariable" begin

    dist = Normal(0, 1)
    name = "x"
    x = RandomVariable(dist, name)

    @test typeof(x) == RandomVariable
    @test x.dist == dist
    @test x.name == name

    @test size(rand(x)) == (1,)
    @test size(rand(x, 10)) == (10,)
end
