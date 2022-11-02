@testset "ResponseSurface" begin
    x = RandomVariable.(Uniform(-π, π), [:x1, :x2, :x3])
    a = 7
    b = 0.1

    ishigami = Model(
        df -> sin.(df.x1) .+ a .* sin.(df.x2) .^ 2 .+ b .* (df.x3 .^ 4) .* sin.(df.x1), :y
    )

    data = sample(x, SobolSampling(100))
    evaluate!(ishigami, data)

    rs = ResponseSurface(data, :y, 6)

    rs_data = select(data, Not(:y))
    evaluate!(rs, rs_data)

    mse = mean((data.y .- rs_data.y) .^ 2)

    @test isapprox(mse, 0; atol=0.02064568539852385)

    @test_throws ErrorException("Degree(p) of ResponseSurface must be non-negative.") ResponseSurface(
        data, :y, -3
    )
end
