using UncertaintyQuantification

# Die normale B-Funktion

inputs = Vector{RandomVariable}(undef, 20)

X = [1,1.1,0.9,1.2,0.8,0.7,1.3,1.4,0.6,0.95]

for i in 1:10
    inputs[i] = RandomVariable(Normal(0, X[i]), Symbol("A$i"))
    inputs[i+10] = RandomVariable(Normal(0, X[i]), Symbol("B$i"))
end


bfkt = Model(df -> df.A1 .* df.A6 + df.A2 .* df.A7 + df.A3 .* df.A8
+ df.A4 .* df.A9 + df.A5 .* df.A10, :fkt)


mc = MonteCarlo(1000)

si = soindex([bfkt], inputs, mc)

#Polyspline der B-Funktion

inputsb = Vector{RandomVariable}(undef, 10)

for i in 1:10
    inputsb[i] = RandomVariable(Normal(0, X[i]), Symbol("A$i"))
end

ps_samples = sample(inputsb, 1000)
evaluate!(bfkt, ps_samples)

ps_samples_ci = convert(Matrix, ps_samples[:,1:10])
ps_samples_f = convert(Matrix, ps_samples[:,11:11])

k = 2

ps_order = df -> [df.A1 df.A2 df.A3 df.A4 df.A5 df.A6 df.A7 df.A8 df.A9 df.A10]

ps_b = createpolyspline(ps_samples_ci, ps_samples_f, k, :fkt, ps_order)

si_ps = soindex([ps_b], inputs, mc)
