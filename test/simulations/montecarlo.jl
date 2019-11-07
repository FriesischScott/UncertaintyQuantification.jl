@testset "MonteCarlo" begin
    mc = MonteCarlo(1000)

    @test isa(mc, MonteCarlo)
    @test mc.n == 1000
end
