@testset "sample" begin

    π = Parameter(3.14, :π)
    x = RandomVariable(Normal(0, 1), :x)
    y = RandomVariable(Normal(0, 1), :y)
    z = RandomVariable(Normal(0, 1), :z)
    jd = JointDistribution([x, y], GaussianCopula([1 0; 0 1]))

    samples = sample([π, jd, z], 10)

    @test size(samples) == (10, 4)
    @test names(samples) == [:π, :x, :y, :z]
end
