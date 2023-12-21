@testset "Distribution Parameters" begin
    # Beta
    d = Beta(0.5, 0.5)
    μ = mean(d)
    σ = std(d)
    @test [distribution_parameters(μ, σ, Beta)...] ≈ [0.5, 0.5]

    # LogNormal
    d = LogNormal(10.0, 5.0)
    μ = mean(d)
    σ = std(d)
    @test [distribution_parameters(μ, σ, LogNormal)...] ≈ [10.0, 5.0]

    # Gumbel
    d = Gumbel(1.0, 2.0)
    μ = mean(d)
    σ = std(d)
    @test [distribution_parameters(μ, σ, Gumbel)...] ≈ [1.0, 2.0]
end
