using UncertaintyQuantification

x = RandomVariable.(Uniform(-π, π), [:x1, :x2, :x3])
a = Parameter(7, :a)
b = Parameter(0.1, :b)

inputs = [x; a; b]

ishigami = Model(df -> sin.(df.x1) .+ df.a .* sin.(df.x2).^2 .+ df.b .* (df.x3.^4) .* sin.(df.x1), :y)

data = sample(inputs, SobolSampling(500))
evaluate!(ishigami, data)

phs = PolyharmonicSpline(data, 2, :y)

si = sobolindices([ishigami], inputs, :y, MonteCarlo(1000))
si_phs = sobolindices([phs], inputs, :y, SobolSampling(500))

println("Sobol indices of the ishigami function: $si")
println("Sobol indices of the polyharmonic spline meta mode: $si_phs")