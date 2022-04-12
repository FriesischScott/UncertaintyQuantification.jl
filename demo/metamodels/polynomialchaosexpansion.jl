using UncertaintyQuantification, Plots

x = RandomVariable.(Uniform(-1, 1), [:x1, :x2])

model = Model(df -> begin
    π .* (df.x1 .- 1) .* sin.(π .* df.x1) .* (1 .- df.x2 .^ 2)
end, :y)

p = 8
Ψ = PolynomialChaosBasis(LegendreBasis.([p, p]), p)

# Estimation by least squares
ls_n = 1000
ls = LeastSquares(SobolSampling(ls_n))
pceLS, samples, mse = polynomialchaos(x, model, Ψ, :y, ls)

println("LS Mean: $(mean(pceLS))")
println("LS Variance: $(var(pceLS))")
println("LS Mean Squared Error: $mse")
println("--------------------------")

# Estimation by full quadrature
gq = GaussQuadrature()
pceGQ, samples = polynomialchaos(x, model, Ψ, :y, gq)

println("GQ Mean: $(mean(pceGQ))")
println("GQ Variance: $(var(pceGQ))")
println("--------------------------")

# Estimation by MC
sampsMC = sample(x, 10^6)
evaluate!(model, sampsMC)

println("MC Mean: $(mean(sampsMC[!,:y]))")
println("MC Variance: $(var(sampsMC[!,:y]))")


# Plots 3 histograms
sampsLS = sample(pceLS, 10^5)
sampsGQ = sample(pceGQ, 10^5)

histogram(sampsLS[!,:y], alpha = 0.5, label = "Least squares", normalize = :probability, bins = 100)
histogram!(sampsGQ[!,:y], alpha = 0.5, label = "Quadrature", normalize = :probability, bins = 100)
histogram!(sampsMC[!,:y], alpha = 0.5, label = "Monte Carlo", normalize = :probability, bins = 100)
