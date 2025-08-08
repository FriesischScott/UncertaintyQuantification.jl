@testset "GaussianMixtureModel" begin
    # Test fitting with a DataFrame and initialization
    df = DataFrame(x1=randn(100), x2=randn(100))
    gmm = GaussianMixtureModel(df, 3; maximum_iterations=1, tolerance=1e-3)

    @test gmm isa RandomUQInput
    @test length(gmm.names) == 2
    @test length(gmm.mixture.components) == 3
    @test dimensions(gmm) == 2
    @test all(names(gmm) .== [:x1, :x2])

    # Test sample generation
    samples = sample(gmm, 10)
    @test nrow(samples) == 10
    @test all(names(samples) .== ["x1", "x2"])

    # Test properties after fitting
    @test all(mean(gmm) .== mean(gmm.mixture))
    @test all(var(gmm) .== var(gmm.mixture))
    @test all(minimum(gmm) == minimum(gmm.mixture))
    @test all(maximum(gmm) == maximum(gmm.mixture))
    @test insupport(gmm, [0.0, 0.0]) == insupport(gmm.mixture, [0.0, 0.0])
    @test pdf(gmm, [0.0, 0.0]) == pdf(gmm.mixture, [0.0, 0.0])
    @test logpdf(gmm, [0.0, 0.0]) == logpdf(gmm.mixture, [0.0, 0.0])

    # Test error handling
    @test_throws ErrorException to_standard_normal_space!(gmm, samples)
    @test_throws ArgumentError fit_gaussian_mixture(3, randn(100, 1))  # Mismatched dimensions
    @test_throws ArgumentError GaussianMixtureModel(MixtureModel([Normal(0, 1)]), [:x1, :x2])  # Mismatched dimensions and names
end
