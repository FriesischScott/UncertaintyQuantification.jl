@testset "Model" begin
    input = DataFrame(a = 1, b = 2)

    model = Model(df -> df.a + 2 * df.b, "c")
    @test typeof(model) == Model

    @test evaluate(model, input) == [5]
    @test model(input) == [5]
end
