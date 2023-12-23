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

            @test mean(ed) ≈ mean(dist) rtol = 0.05
            @test var(ed) ≈ var(dist) rtol = 0.05
            @test minimum(ed) == minimum(ed.data)
            @test maximum(ed) == maximum(ed.data)

            samples = range(minimum(ed), maximum(ed); length=100)

            @test all(insupport.(ed, samples))
            @test cdf.(ed, samples) ≈ cdf.(dist, samples) atol = 0.05

            @test all(pdf.(ed, samples) .>= 0)
            pdf_area, _ = hquadrature(x -> pdf(ed, x), minimum(ed), maximum(ed))

            @test pdf_area ≈ 1 atol = 1e-3
            @test logpdf.(ed, samples) ≈ log.(pdf.(ed, samples))

            r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
            @test quantile.(ed, r) ≈ quantile.(dist, r) atol = 0.05
        end
    end
end
