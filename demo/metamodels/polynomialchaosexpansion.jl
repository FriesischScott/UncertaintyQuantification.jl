using UncertaintyQuantification, DataFrames

x = RandomVariable.(Uniform(-1, 1), [:x1, :x2])

model = Model(df -> π .* (df.x1 .- 1) .* sin.(π .* df.x1) .* (1 .- df.x2 .^ 2), :y)

data = sample(x, SobolSampling(100))
evaluate!(model, data)

pce = PolynomialChaosExpansion(data, x, 4, :y)

test_data = select(data, Not(:y))
evaluate!(pce, test_data)

abs.(data.y .- test_data.y)
