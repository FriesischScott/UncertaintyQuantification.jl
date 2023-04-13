@testset "FORM" begin
    form = FORM()

    x1 = RandomVariable(Normal(200, 20), :x1)
    x2 = RandomVariable(Normal(150, 10), :x2)

    inputs = [x1, x2]

    model = Model(df -> df.x1 .- df.x2, :y)

    pf, β, dp = probability_of_failure(df -> df.x1 .- df.x2, inputs, form)

    @test round(pf; digits=4) ≈ 0.0127
    @test round(β; digits=4) ≈ 2.2361

    pf, β, dp = probability_of_failure(model, df -> df.y, inputs, form)

    @test round(pf; digits=4) ≈ 0.0127
    @test round(β; digits=4) ≈ 2.2361
end
