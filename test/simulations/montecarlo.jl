@testset "MonteCarlo" begin
    mc = MonteCarlo(1000)

    @test typeof(mc) == MonteCarlo
    @test mc.n == 1000
end
