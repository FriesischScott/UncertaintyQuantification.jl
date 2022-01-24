using UncertaintyQuantification, DataFrames

x = RandomVariable.(Uniform(-π, π), [:x1, :x2, :x3])
a = Parameter(7, :a)
b = Parameter(0.05, :b)

inputs = [x; a; b]

ishigami = Model(df -> sin.(df.x1) .+ df.a .* sin.(df.x2).^2 .+ df.b .* (df.x3.^4) .* sin.(df.x1), :y)

data = sample(inputs, MonteCarlo(1000))
evaluate!(ishigami, data)

pce = PolynomialChaosExpansion(data, inputs, 3, :y)

test_data = select(data, Not(:y))
evaluate!(pce, test_data)

abs.(data.y .- test_data.y)