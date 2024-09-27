@testset "Spectral Representation" begin
    N = 128
    ω_u = 4π
    Δω = ω_u / N
    T0 = 2π/Δω
    Δt = 2π/(2*ω_u)

    t = collect(Δt:Δt:(T0))
    ω = collect(0:Δω:(ω_u - Δω))

    sd = ShinozukaDeodatis(ω, 1, 1)
    Sval = evaluate(sd)

    srm_obj = SpectralRepresentation(sd, t, :ShnzkSR)
    ϕ = sample(srm_obj)
    ϕ_temp = copy(ϕ)
    x = evaluate(srm_obj, ϕ[1, :])

    Ŝ = Δt^2 / T0 * abs.(dft(x)) .^ 2 / (2π)
    S̃ = periodogram(x,t)
    S̄ = periodogram(x,t,false)

    x_2 = srm_obj(collect(ϕ[1, :]))
    x_3 = srm_obj(ϕ[1, :])

    to_standard_normal_space!(srm_obj, ϕ)
    to_physical_space!(srm_obj, ϕ)

    @test Sval ≈ Ŝ[1:N]
    @test Sval ≈ S̃[1:N]
    @test Ŝ ≈ S̃
    @test S̄ ≈ 2 .* Ŝ[1:N]
    @test S̄ ≈ 2 .* S̃[1:N]  

    @test x ≈ x_2
    @test x ≈ x_3
    @test x_2 ≈ x_3

    @test ϕ ≈ ϕ_temp

end