using UncertaintyQuantification

σx = [1, 1.1, 0.9, 1.2, 0.8]
σω = [0.7, 1.3, 1.4, 0.6, 0.95]

X = RandomVariable.(Normal.(0, σx), [:X1, :X2, :X3, :X4, :X5])
ω = RandomVariable.(Normal.(0, σω), [:ω1, :ω2, :ω3, :ω4, :ω5])

B = Model(df -> df.X1 .* df.ω1
    + df.X2 .* df.ω2
    + df.X3 .* df.ω3
    + df.X4 .* df.ω4
    + df.X5 .* df.ω5, :B)

mc = MonteCarlo(100000)
sobol = SobolSampling(100000)

si_mc = sobolindices([B], [X; ω], :B, mc)
si_sobol = sobolindices([B], [X; ω], :B, sobol)
