@testset "SpectralRepresentation" begin
    ω_u = 4 * π
    N = 128
    Δω = ω_u / N
    T0 = 64
    dt = 0.25

    t = collect(0:dt:(T0 - dt))
    ω = collect(0:Δω:(ω_u - Δω))

    sd = ShinozukaDeodatis(ω, 1, 1)
    Sval = evaluate(sd)

    srm_obj = SpectralRepresentation(sd, t, :ShnzkSR)
    ϕ = sample(srm_obj)
    x = evaluate(srm_obj, ϕ[1, :])

    S_hat = dt^2 / T0 * abs.(dft(x)) .^ 2 / (2 * π)

    @test Sval≈S_hat[1:N]
end