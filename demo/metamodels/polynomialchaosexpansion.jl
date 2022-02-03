using UncertaintyQuantification, DataFrames

x = RandomVariable.(Uniform(-1, 1), [:x1, :x2])

model = Model(df -> begin
    π .* (df.x1 .- 1) .* sin.(π .* df.x1) .* (1 .- df.x2 .^ 2)
end, :y)

data = sample(x, SobolSampling(1000))
evaluate!(model, data)

Ψ = LegendreBasis(10, 2)
pce = PolynomialChaosExpansion(data, x, Ψ, :y)

μ = pce.y[1]
global σ² = 0.0

for i in 2:length(pce.y)
    global σ² +=
        (pce.y[i] * sqrt(1 / (2 * pce.Ψ.indices[i][1] + 1)) * sqrt(1 / (2 * pce.Ψ.indices[i][2] + 1)))^2
end

println("Mean: $μ")
println("Variance: $σ²")

test_data = select(data, Not(:y))
evaluate!(pce, test_data)

ϵ = sum(abs.(data.y .- test_data.y))
println("Residual sum: $ϵ")
