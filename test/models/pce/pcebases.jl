@testset "PolyChaosBasis" begin
    @testset "LegendreBasis" begin
        x = range(-1, 1; length=10)

        Ψ = LegendreBasis(4, false)

        @test evaluate(Ψ, x, 0) == ones(10)
        @test evaluate(Ψ, x, 1) == x
        @test evaluate(Ψ, x, 2) ≈ (1 / 2) * (3x .^ 2 .- 1)
        @test evaluate(Ψ, x, 3) ≈ (1 / 2) * (5x .^ 3 .- 3x)
        @test evaluate(Ψ, x, 4) ≈ (1 / 8) * (35x .^ 4 .- 30x .^ 2 .+ 3)

        Ψ = LegendreBasis(4, true)

        @test evaluate(Ψ, x, 0) == ones(10)
        @test evaluate(Ψ, x, 1) == x .* sqrt(3)
        @test evaluate(Ψ, x, 2) ≈ (1 / 2) * (3x .^ 2 .- 1) .* sqrt(5)
        @test evaluate(Ψ, x, 3) ≈ (1 / 2) * (5x .^ 3 .- 3x) .* sqrt(7)
        @test evaluate(Ψ, x, 4) ≈ (1 / 8) * (35x .^ 4 .- 30x .^ 2 .+ 3) .* sqrt(9)

        x, w = gausslegendre(10)
        @test quadrature_nodes(10, Ψ) == x
        @test quadrature_weights(10, Ψ) == w ./ 2
    end

    @testset "HermiteBasis" begin
        x = range(-3, 3; length=10)

        Ψ = HermiteBasis(4, false)

        @test evaluate(Ψ, x, 0) == ones(10)
        @test evaluate(Ψ, x, 1) == x
        @test evaluate(Ψ, x, 2) ≈ x .^ 2 .- 1
        @test evaluate(Ψ, x, 3) ≈ x .^ 3 .- 3x
        @test evaluate(Ψ, x, 4) ≈ x .^ 4 .- 6x .^ 2 .+ 3

        Ψ = HermiteBasis(4, true)

        @test evaluate(Ψ, x, 0) == ones(10)
        @test evaluate(Ψ, x, 1) == x
        @test evaluate(Ψ, x, 2) ≈ (x .^ 2 .- 1) / sqrt(2)
        @test evaluate(Ψ, x, 3) ≈ (x .^ 3 .- 3x) / sqrt(6)
        @test evaluate(Ψ, x, 4) ≈ (x .^ 4 .- 6x .^ 2 .+ 3) / sqrt(24)

        x, w = gausshermite(10)
        @test quadrature_nodes(10, Ψ) == sqrt(2) .* x
        @test quadrature_weights(10, Ψ) == w ./ sqrt(π)
    end

    @testset "multivariate_indices" begin
        @test multivariate_indices(2, 1) == [[0], [1], [2]]
        @test multivariate_indices(2, 2) == [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]]
        @test multivariate_indices(2, 3) == [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [2, 0, 0],
            [1, 1, 0],
            [1, 0, 1],
            [0, 2, 0],
            [0, 1, 1],
            [0, 0, 2],
        ]
    end
end
