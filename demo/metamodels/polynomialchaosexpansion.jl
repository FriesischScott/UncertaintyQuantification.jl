using UncertaintyQuantification

x = RandomVariable.(Uniform(-1, 1), [:x1, :x2])

model = Model(df -> begin
    π .* (df.x1 .- 1) .* sin.(π .* df.x1) .* (1 .- df.x2 .^ 2)
end, :y)

Ψ = LegendreBasis(8, 2)

# Estimation by least squares
ls = LeastSquares(SobolSampling(1000))
pce, samples, mse = polynomialchaos(x, [model], Ψ, :y, ls)

println("LS Mean: $(mean(pce))")
println("LS Variance: $(var(pce))")
println("LS Mean Squared Error: $mse")

# Estimation by full quadrature
gq = GaussQuadrature()
pce, samples = polynomialchaos(x, [model], Ψ, :y, gq)

println("GQ Mean: $(mean(pce))")
println("GQ Variance: $(var(pce))")

# Estimation by sparse quadrature
sq = SparseQuadrature()
pce, samples = polynomialchaos(x, [model], Ψ, :y, sq)

println("SQ Mean: $(mean(pce))")
println("SQ Variance: $(var(pce))")
