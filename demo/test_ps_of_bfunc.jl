using UncertaintyQuantification

# Die normale B-Funktion

σx = [1, 1.1, 0.9, 1.2, 0.8]
σω = [0.7, 1.3, 1.4, 0.6, 0.95]

X = RandomVariable.(Normal.(0, σx), [:X1, :X2, :X3, :X4, :X5])
ω = RandomVariable.(Normal.(0, σω), [:ω1, :ω2, :ω3, :ω4, :ω5])

B = Model(df -> df.X1 .* df.ω1 
    + df.X2 .* df.ω2
    + df.X3 .* df.ω3
    + df.X4 .* df.ω4
    + df.X5 .* df.ω5, :B)

mc = MonteCarlo(10000)

si = sobolindices([B], [X; ω], :B, mc)

# Polyspline der B-Funktion

inputsb = Vector{RandomVariable}(undef, 10)

for i in 1:5
    inputsb[i] = RandomVariable(Normal(0, X[i]), Symbol("x$i"))
    inputsb[i + 5] = RandomVariable(Normal(0, X[i + 5]), Symbol("w$i"))
end

ps_samples = sample(inputsb, 1000)
evaluate!(bfkt, ps_samples)

ps_samples_ci = convert(Matrix, ps_samples[:,1:10])
ps_samples_f = convert(Matrix, ps_samples[:,11:11])

k = 2

ps_order = df -> [df.x1 df.x2 df.x3 df.x4 df.x5 df.w1 df.w2 df.w3 df.w4 df.w5]

ps_b = PolyharmonicSpline(ps_samples_ci, ps_samples_f, k, :fkt)

si_ps_k_2 = sobolindices([ps_b], inputs, mc)
