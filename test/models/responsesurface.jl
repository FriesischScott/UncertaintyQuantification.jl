using UncertaintyQuantification 
using Test
using DataFrames

@testset "ResponseSurface" begin
    x = RandomVariable.(Uniform(-π ,π), [:x1, :x2, :x3])
    a = 7
    b = 0.1

    ishigami = Model(
        df ->
            sin.(df.x1) .+ a .* sin.(df.x2) .^ 2 .+ b .* (df.x3 .^ 4) .* sin.(df.x1),
        :y,
    )

    data = sample(x, 100)

    evaluate!(ishigami, data)

    rs = ResponseSurface(data, :y, 6)

    rs_data = select(data, Not(:y))
    evaluate!(rs, rs_data)

    mse = mean((data.y .- rs_data.y) .^ 2)

    @test isapprox(mse, 0; atol=eps(Float64))
end