using UncertaintyQuantification, DataFrames

x = RandomVariable.(Uniform(-1, 1), [:x1, :x2])

model = Model(df -> π .* (df.x1 .- 1) .* sin.(π .* df.x1) .* (1 .- df.x2 .^ 2), :y)

data = sample(x, SobolSampling(10000))
evaluate!(model, data)

pce = PolynomialChaosExpansion(data, x, 4, :y)

μ = pce.y[1]
σ² = sum(pce.y[2:end] .^ 2)

println("Mean: $μ")
println("Variance: $σ²")

test_data = select(data, Not(:y))
evaluate!(pce, test_data)

ϵ = sum(abs.(data.y .- test_data.y))
println("Residual sum: $ϵ")
