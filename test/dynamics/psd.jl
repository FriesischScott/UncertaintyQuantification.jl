@testset "PSD generation" begin
    ω = collect(0:0.6:150)
    cp = CloughPenzien(ω, 0.1, 0.8π, 0.6, 8π, 0.6)
    kt = KanaiTajimi(ω, 0.25, 5, 0.75)
    sd = ShinozukaDeodatis(ω, 1, 1)
    ep = EmpiricalPSD(ω, evaluate(cp))

    @test isa(cp, CloughPenzien)
    @test isa(kt, KanaiTajimi)
    @test isa(sd, ShinozukaDeodatis)
    @test evaluate(cp) ≈ evaluate(ep)

end