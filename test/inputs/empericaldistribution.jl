@testset "EmpiricalDistribution" begin
    n = 10^5

    distributions = [Normal(1, 1), Rayleigh(1), Gumbel()]
    names = ["Normal", "Rayleigh", "Gumbel"]
    for (dist, name) in zip(distributions, names)
        @testset "$name" begin
            x = rand(dist, n)
            ed = EmpiricalDistribution(x)
            y = rand(ed, n)

            h0 = ApproximateTwoSampleKSTest(x, y)
            @test pvalue(h0) > 0.05

            @test ≈(mean(ed), mean(dist); rtol=0.05)
            @test ≈(var(ed), var(dist); rtol=0.05)
            @test minimum(ed) == minimum(ed.data)
            @test maximum(ed) == maximum(ed.data)

            samples = rand(Uniform(minimum(ed), maximum(ed)), 100)
            @test all(insupport.(ed, samples))
            @test all(isapprox.(cdf.(ed, samples), cdf.(dist, samples); atol=0.05))
            @test all(isapprox.(pdf.(ed, samples), pdf.(dist, samples); atol=0.05))
            @test all(isapprox.(logpdf.(ed, samples), log.(pdf.(ed, samples)); atol=0.05))

            r = [0.25, 0.5, 0.75]
            @test all(isapprox.(quantile.(ed, r), quantile.(dist, r); atol=0.05))
        end
    end
end
