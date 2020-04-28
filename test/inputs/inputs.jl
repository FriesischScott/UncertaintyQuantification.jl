π = Parameter(3.14, :π)
x = RandomVariable(Normal(0, 1), :x)
y = RandomVariable(Normal(1, 1), :y)
z = RandomVariable(Normal(0, 1), :z)
jd = JointDistribution([x, y], GaussianCopula([1 0; 0 1]))

inputs = [π, jd, z]

@testset "Inputs" begin
    @testset "sample" begin
        samples = sample(inputs, 10)

        @test size(samples) == (10, 4)
        @test names(samples) == [:π, :x, :y, :z]
    end

    @testset "mean" begin
        means = mean(inputs)
        @test size(means) == (1, 4)
        @test means.π[1] == 3.14
        @test means.x[1] == 0
        @test means.y[1] == 1
        @test means.z[1] == 0
    end

    @testset "names" begin
        @test names(inputs) == [:π, :x, :y, :z]
    end
end