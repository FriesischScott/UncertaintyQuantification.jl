@testset "EmpiricalDistribution" begin
    data = [
        21.370,
        19.435,
        20.363,
        20.632,
        20.404,
        19.893,
        21.511,
        19.905,
        22.018,
        19.93,
        31.304,
        32.286,
        28.611,
        29.721,
        29.866,
        30.635,
        29.715,
        27.343,
        27.559,
        31.32,
        39.693,
        38.218,
        39.828,
        41.214,
        41.895,
        39.569,
        39.742,
        38.236,
        40.460,
        39.36,
        50.455,
        50.704,
        51.035,
        49.391,
        50.504,
        48.282,
        49.215,
        49.149,
        47.585,
        50.03,
    ]

    ed = EmpiricalDistribution(data)

    @test ed.h ≈ 1.0870436775976158

    @test mean(ed) ≈ mean(data)
    @test var(ed) ≈ var(data)

    samples = range(15, 55; length=250)

    @test all(insupport.(ed, samples))
    @test all(pdf.(ed, samples) .>= 0)

    pdf_area, _ = hquadrature(x -> pdf(ed, x), 15, 55)

    @test pdf_area ≈ 1 atol = 0.01
    @test logpdf.(ed, samples) ≈ log.(pdf.(ed, samples))
end
