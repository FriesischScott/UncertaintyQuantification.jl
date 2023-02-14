@testset "EmpiricalDistribution" begin
    n = 10^4

    distributions = [Normal(), Exponential()]
    names = ["Normal", "Exponential"]

    for (dist, name) in zip(distributions, names)
        @testset "$name" begin
            x = rand(dist, n)

            ed = EmpiricalDistribution(x)
            y = rand(ed, n)

            h0 = ApproximateTwoSampleKSTest(x, y)
            @test pvalue(h0) > 0.05
        end
    end
end
